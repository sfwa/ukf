/*
Reading a recording: the ".pack" files written by explainer/tools/pack_euroc.py.
That file's header comment is the description of the format.

The IMU values come back as the same doubles which the benchmark reads from
the dataset's CSV files, bit for bit (the tool that writes the packs proves
that for every value).
*/

const MAGIC = 0x58464b55;   /* "UKFX", little-endian */
const FORMAT_VERSION = 1;
const ITEM_SIZE = {uint8: 1, int16: 2, int32: 4};

async function gunzip(buffer) {
    const stream = new Blob([buffer]).stream().pipeThrough(new DecompressionStream('gzip'));
    return new Response(stream).arrayBuffer();
}

/*
One array as plain numbers, "components" per row. Undoes the three storage
steps in reverse: bytes by position, values by component, then differences.
*/
function decodeIntegers(bytes, description) {
    const {type, count, components, delta} = description;
    const size = ITEM_SIZE[type];
    const total = count * components;
    const start = description.offset;
    const signed = type !== 'uint8';
    const limit = 2 ** (8 * size);
    const out = new Float64Array(total);
    for (let component = 0; component < components; component++) {
        for (let row = 0; row < count; row++) {
            const item = component * count + row;
            let value = 0;
            for (let byte = size - 1; byte >= 0; byte--) {
                value = value * 256 + bytes[start + byte * total + item];
            }
            if (signed && value >= limit / 2) {
                value -= limit;
            }
            out[row * components + component] = value;
        }
    }
    for (let pass = 0; pass < delta; pass++) {
        for (let i = components; i < total; i++) {
            out[i] += out[i - components];
        }
    }
    return out;
}

/* Spread an array kept at every Nth sample over all samples, linearly. */
function spread(values, components, every, samples) {
    if (every === 1) {
        return values;
    }
    const rows = values.length / components;
    const out = new Float64Array(samples * components);
    for (let i = 0; i < samples; i++) {
        const position = i / every;
        const left = Math.min(Math.floor(position), rows - 1);
        const right = Math.min(left + 1, rows - 1);
        const mix = position - left;
        for (let k = 0; k < components; k++) {
            out[i * components + k] = values[left * components + k] * (1 - mix) + values[right * components + k] * mix;
        }
    }
    return out;
}

/*
Decode a pack file into a recording:

    name, samples, duration, header
    time            seconds since the first sample
    gyroscope       rad/s, IMU frame, 3 per sample
    accelerometer   m/s^2, IMU frame, 3 per sample
    reference       valid (0 or 1), attitude (w x y z), position, velocity,
                    gyroBias, accelerometerBias: all per IMU sample
*/
export async function decodePack(buffer) {
    let bytes = new Uint8Array(buffer);
    /* A server may already have removed the gzip layer. */
    if (bytes[0] === 0x1f && bytes[1] === 0x8b) {
        bytes = new Uint8Array(await gunzip(buffer));
    }
    const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
    if (view.getUint32(0, true) !== MAGIC || view.getUint32(4, true) !== FORMAT_VERSION) {
        throw new Error('This is not a recording which this page can read');
    }
    const headerLength = view.getUint32(8, true);
    const header = JSON.parse(new TextDecoder().decode(bytes.subarray(12, 12 + headerLength)));
    const body = bytes.subarray(12 + headerLength);
    const samples = header.samples;
    const read = name => decodeIntegers(body, header.arrays[name]);
    const scaled = name => {
        const description = header.arrays[name];
        const values = read(name);
        for (let i = 0; i < values.length; i++) {
            values[i] *= description.scale;
        }
        return spread(values, description.components, description.every, samples);
    };

    /* The decode rules of the pack tool. Their order of operations matters. */
    const time = read('time_ns');
    for (let i = 0; i < samples; i++) {
        time[i] *= 1e-9;
    }
    const gyroscope = read('gyroscope');
    for (let i = 0; i < gyroscope.length; i++) {
        gyroscope[i] = gyroscope[i] / 25 / 180 * Math.PI;
    }
    const accelerometer = read('accelerometer');
    for (let i = 0; i < accelerometer.length; i++) {
        accelerometer[i] = accelerometer[i] / 1200 * 9.80665;
    }

    /* Rounding leaves the reference attitude very slightly off unit length. */
    const attitude = scaled('reference_attitude');
    for (let i = 0; i < samples; i++) {
        const norm = Math.hypot(attitude[4 * i], attitude[4 * i + 1], attitude[4 * i + 2], attitude[4 * i + 3]);
        for (let k = 0; k < 4; k++) {
            attitude[4 * i + k] /= norm;
        }
    }

    return {
        name: header.sequence,
        samples,
        duration: time[samples - 1],
        header,
        time,
        gyroscope,
        accelerometer,
        reference: {
            valid: Uint8Array.from(read('reference_valid')),
            attitude,
            position: scaled('reference_position'),
            velocity: scaled('reference_velocity'),
            gyroBias: scaled('reference_gyro_bias'),
            accelerometerBias: scaled('reference_accelerometer_bias'),
        },
    };
}

export async function fetchRecording(url) {
    const response = await fetch(url);
    if (!response.ok) {
        throw new Error(`Cannot load ${url} (${response.status})`);
    }
    return decodePack(await response.arrayBuffer());
}
