/*
The aircraft drawn in the three-dimensional figures: a six-rotor
multicopter with a camera bar at the front, like the AscTec Firefly which
flew the EuRoC recordings. The shape is a plain sketch, not a scale model.

The model is built in the vehicle frame: X forward, Y left, Z up along the
rotor thrust. The IMU sits in the camera bar and is tilted, so its own
frame is different: thrust is 20 degrees from the IMU's +X axis. The
estimator's attitude is that of the IMU frame, and imuToVehicle() gives
the fixed turn between the two from the estimator's own constants.
*/

import {box, cylinder, lathe, sphere, arrowShape, transform, merge} from './gl.js';
import {mMul, mTranslate, mFromQuat, mFromAxes, mScale, mPose, qFromAxisAngle, qBetween, sub, norm, scale, DEG} from './math.js';

const ROTOR_RADIUS = 0.215;     /* centre of the vehicle to each rotor hub (m) */
const PROPELLER_RADIUS = 0.1;
export const VEHICLE_FOOT = 0.13;      /* the feet are this far below the vehicle's origin (m) */

/* A thin cylinder from one point to another. */
function strut(from, to, radius) {
    const direction = sub(to, from);
    const shape = cylinder(radius, norm(direction), 8);
    return transform(shape, mMul(mTranslate(from), mFromQuat(qBetween([0, 0, 1], direction))));
}

function buildShapes() {
    const frame = [], motors = [], discs = [], blades = [];

    /* Centre: a rounded hub with a box for the battery under it. */
    frame.push(lathe([[0, -0.034], [0.068, -0.034], [0.08, -0.012], [0.08, 0.012], [0.056, 0.03], [0, 0.034]], 28));
    frame.push(transform(box(0.15, 0.064, 0.036), mTranslate([-0.01, 0, -0.052])));

    for (let k = 0; k < 6; k++) {
        const angle = (30 + 60 * k) * DEG;
        const turn = mFromQuat(qFromAxisAngle([0, 0, 1], angle));
        const hub = [ROTOR_RADIUS * Math.cos(angle), ROTOR_RADIUS * Math.sin(angle), 0];
        frame.push(transform(box(ROTOR_RADIUS, 0.018, 0.012), mMul(turn, mTranslate([ROTOR_RADIUS / 2, 0, 0]))));
        motors.push(transform(cylinder(0.021, 0.03, 16), mTranslate([hub[0], hub[1], 0.004])));
        motors.push(transform(cylinder(0.007, 0.014, 10), mTranslate([hub[0], hub[1], 0.034])));
        discs.push(transform(cylinder(PROPELLER_RADIUS, 0.002, 36), mTranslate([hub[0], hub[1], 0.042])));
        /* A stopped propeller: two blades, each at its own angle so they do not look alike. */
        const bladeTurn = mFromQuat(qFromAxisAngle([0, 0, 1], angle + (35 + 50 * k) * DEG));
        blades.push(transform(box(2 * PROPELLER_RADIUS, 0.02, 0.003), mMul(mTranslate([hub[0], hub[1], 0.043]), bladeTurn)));
    }

    /* Four legs with round feet. */
    for (const [sx, sy] of [[1, 1], [1, -1], [-1, 1], [-1, -1]]) {
        const top = [0.05 * sx, 0.045 * sy, -0.03], foot = [0.105 * sx, 0.115 * sy, -VEHICLE_FOOT + 0.012];
        frame.push(strut(top, foot, 0.005));
        frame.push(transform(sphere(0.012, 10), mTranslate(foot)));
    }

    /* The camera bar at the front, which also holds the IMU. It looks 20 degrees down. */
    const pitch = mFromQuat(qFromAxisAngle([0, 1, 0], 20 * DEG));
    const mount = mMul(mTranslate([0.108, 0, -0.04]), pitch);
    const sensor = [transform(box(0.028, 0.15, 0.038), mount)];
    for (const side of [-1, 1]) {
        const lens = transform(cylinder(0.013, 0.014, 14), mFromQuat(qFromAxisAngle([0, 1, 0], 90 * DEG)));
        sensor.push(transform(lens, mMul(mount, mTranslate([0.012, 0.055 * side, 0]))));
    }
    frame.push(strut([0.06, 0, -0.03], [0.1, 0, -0.04], 0.008));

    return {frame: merge(frame), motors: merge(motors), discs: merge(discs), blades: merge(blades), sensor: merge(sensor)};
}

/*
The fixed turn from the vehicle frame to the IMU frame, as a matrix: its
columns are the vehicle's forward, left and up directions in IMU
coordinates. "constants" are the estimator module's constants.
*/
export function vehicleInImu(constants) {
    const up = [constants.thrust_axis_x, constants.thrust_axis_y, constants.thrust_axis_z];
    const forward = [constants.lateral_axis_2_x, constants.lateral_axis_2_y, constants.lateral_axis_2_z];
    const left = scale([constants.lateral_axis_1_x, constants.lateral_axis_1_y, constants.lateral_axis_1_z], -1);
    return mFromAxes(forward, left, up);
}

export class Vehicle {
    constructor(renderer, constants) {
        const shapes = buildShapes();
        this.meshes = Object.fromEntries(Object.entries(shapes).map(([name, shape]) => [name, renderer.mesh(shape)]));
        this.mount = constants ? vehicleInImu(constants) : mFromAxes([1, 0, 0], [0, 1, 0], [0, 0, 1]);
    }

    /*
    Draw the vehicle. "attitude" turns the IMU frame into the world frame,
    as the estimator reports it. Options: size, spinning (rotors turning),
    tint and opacity (to draw a see-through copy in one colour).
    */
    draw(renderer, palette, position, attitude, options = {}) {
        const {size = 1, spinning = false, tint = null, opacity = 1} = options;
        const model = mMul(mPose(position, attitude, size), this.mount);
        const color = name => tint ?? palette[name];
        renderer.draw(this.meshes.frame, model, color('bodyRgb'), {opacity});
        renderer.draw(this.meshes.motors, model, color('body2Rgb'), {opacity});
        renderer.draw(this.meshes.sensor, model, color('body2Rgb'), {opacity});
        if (spinning) {
            renderer.draw(this.meshes.discs, model, color('body3Rgb'), {opacity: 0.3 * opacity});
        } else {
            renderer.draw(this.meshes.blades, model, color('body3Rgb'), {opacity});
        }
    }
}

let arrowMesh = null;
let ballMesh = null;

/* An arrow from a world point along a vector (its length is the vector's). */
export function drawArrow(renderer, from, vector, color, options = {}) {
    const length = norm(vector);
    if (length < 1e-6) {
        return;
    }
    arrowMesh ??= renderer.mesh(arrowShape());
    const thickness = options.thickness ?? 1;
    /* Keep the arrow's width the same whatever its length. */
    const model = mMul(mMul(mTranslate(from), mFromQuat(qBetween([0, 0, 1], vector))), mScale([thickness, thickness, length]));
    renderer.draw(arrowMesh, model, color, {opacity: options.opacity ?? 1});
}

/* A small ball at a world point. */
export function drawBall(renderer, at, radius, color, options = {}) {
    ballMesh ??= renderer.mesh(sphere(1, 16));
    renderer.draw(ballMesh, mMul(mTranslate(at), mScale(radius)), color, options);
}
