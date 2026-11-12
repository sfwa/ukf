/*
A small WebGL renderer for the three-dimensional figures: shaded solid
shapes only.

A browser allows a page only a few WebGL canvases, and this page has many
figures. So all of them share one hidden WebGL canvas. A figure draws its
solids into it and then copies the picture onto its own ordinary canvas,
where it can also draw lines and labels under and over them.
*/

const VERTEX_SHADER = `
attribute vec3 position;
attribute vec3 normal;
uniform mat4 viewProjection;
uniform mat4 model;
varying vec3 worldNormal;
void main() {
    worldNormal = mat3(model) * normal;
    gl_Position = viewProjection * model * vec4(position, 1.0);
}`;

const FRAGMENT_SHADER = `
precision mediump float;
uniform vec4 color;
uniform vec3 light;
uniform vec3 skyLight;
uniform vec3 groundLight;
uniform float keyLight;
uniform float flatShade;
varying vec3 worldNormal;
void main() {
    vec3 n = normalize(worldNormal);
    if (!gl_FrontFacing) {
        n = -n;
    }
    /* Soft light from the sky and the ground, plus one directional light. */
    vec3 ambient = mix(groundLight, skyLight, n.z * 0.5 + 0.5);
    float direct = max(dot(n, light), 0.0);
    vec3 lit = color.rgb * mix(ambient + keyLight * direct, vec3(1.0), flatShade);
    gl_FragColor = vec4(lit * color.a, color.a);
}`;

let shared = null;

export class Renderer {
    /* The one renderer of the page, or null if WebGL is not available. */
    static get() {
        if (shared === null) {
            try {
                shared = new Renderer();
            } catch (error) {
                console.warn('WebGL is not available; three-dimensional figures are not drawn.', error);
                shared = false;
            }
        }
        return shared || null;
    }

    constructor() {
        this.canvas = document.createElement('canvas');
        const gl = this.canvas.getContext('webgl', {antialias: true, alpha: true, premultipliedAlpha: true});
        if (!gl) {
            throw new Error('No WebGL context');
        }
        this.gl = gl;
        const program = gl.createProgram();
        for (const [type, source] of [[gl.VERTEX_SHADER, VERTEX_SHADER], [gl.FRAGMENT_SHADER, FRAGMENT_SHADER]]) {
            const shader = gl.createShader(type);
            gl.shaderSource(shader, source);
            gl.compileShader(shader);
            if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
                throw new Error(gl.getShaderInfoLog(shader));
            }
            gl.attachShader(program, shader);
        }
        gl.bindAttribLocation(program, 0, 'position');
        gl.bindAttribLocation(program, 1, 'normal');
        gl.linkProgram(program);
        if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
            throw new Error(gl.getProgramInfoLog(program));
        }
        this.program = program;
        this.uniforms = {};
        for (const name of ['viewProjection', 'model', 'color', 'light', 'skyLight', 'groundLight', 'keyLight', 'flatShade']) {
            this.uniforms[name] = gl.getUniformLocation(program, name);
        }
        this.width = 0;
        this.height = 0;
    }

    /* Make a mesh from {positions, normals}: three numbers per corner, three corners per triangle. */
    mesh(shape) {
        const gl = this.gl;
        const data = new Float32Array(shape.positions.length * 2);
        for (let i = 0; i < shape.positions.length / 3; i++) {
            data.set(shape.positions.subarray(3 * i, 3 * i + 3), 6 * i);
            data.set(shape.normals.subarray(3 * i, 3 * i + 3), 6 * i + 3);
        }
        const buffer = gl.createBuffer();
        gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
        gl.bufferData(gl.ARRAY_BUFFER, data, gl.STATIC_DRAW);
        return {buffer, count: shape.positions.length / 3};
    }

    /*
    Start a picture of the given size in CSS pixels. "lighting" may set
    light (a direction towards the light), sky, ground and key.
    */
    begin(width, height, dpr, viewProjection, lighting = {}) {
        const gl = this.gl;
        this.width = Math.round(width * dpr);
        this.height = Math.round(height * dpr);
        /* Grow only, so that figures of different sizes do not keep resizing the canvas. */
        if (this.canvas.width < this.width || this.canvas.height < this.height) {
            this.canvas.width = Math.max(this.canvas.width, this.width);
            this.canvas.height = Math.max(this.canvas.height, this.height);
        }
        gl.viewport(0, 0, this.width, this.height);
        gl.enable(gl.SCISSOR_TEST);
        gl.scissor(0, 0, this.width, this.height);
        gl.clearColor(0, 0, 0, 0);
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
        gl.enable(gl.DEPTH_TEST);
        gl.enable(gl.BLEND);
        gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
        gl.useProgram(this.program);
        gl.uniformMatrix4fv(this.uniforms.viewProjection, false, viewProjection);
        const light = lighting.light ?? [0.35, -0.45, 0.82];
        const length = Math.hypot(...light);
        gl.uniform3f(this.uniforms.light, light[0] / length, light[1] / length, light[2] / length);
        gl.uniform3fv(this.uniforms.skyLight, lighting.sky ?? [0.62, 0.62, 0.64]);
        gl.uniform3fv(this.uniforms.groundLight, lighting.ground ?? [0.34, 0.33, 0.32]);
        gl.uniform1f(this.uniforms.keyLight, lighting.key ?? 0.5);
    }

    /* Draw a mesh. "color" is [r, g, b] from 0 to 1. Options: opacity, flat (no shading). */
    draw(mesh, model, color, options = {}) {
        const gl = this.gl;
        const opacity = options.opacity ?? 1;
        gl.depthMask(opacity >= 1);
        gl.uniformMatrix4fv(this.uniforms.model, false, model);
        gl.uniform4f(this.uniforms.color, color[0], color[1], color[2], opacity);
        gl.uniform1f(this.uniforms.flatShade, options.flat ? 1 : 0);
        gl.bindBuffer(gl.ARRAY_BUFFER, mesh.buffer);
        gl.enableVertexAttribArray(0);
        gl.enableVertexAttribArray(1);
        gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 24, 0);
        gl.vertexAttribPointer(1, 3, gl.FLOAT, false, 24, 12);
        gl.drawArrays(gl.TRIANGLES, 0, mesh.count);
    }

    /* Copy the finished picture onto a figure's canvas (whose transform is in CSS pixels). */
    finish(ctx, width, height) {
        this.gl.depthMask(true);
        ctx.drawImage(this.canvas, 0, this.canvas.height - this.height, this.width, this.height, 0, 0, width, height);
    }
}

/* Shapes. Each returns {positions, normals} as Float32Arrays. */

function build(triangles) {
    const positions = new Float32Array(triangles.length * 9);
    const normals = new Float32Array(triangles.length * 9);
    triangles.forEach((triangle, t) => {
        for (let corner = 0; corner < 3; corner++) {
            positions.set(triangle.points[corner], 9 * t + 3 * corner);
            normals.set(triangle.normals[corner], 9 * t + 3 * corner);
        }
    });
    return {positions, normals};
}

function faceNormal(a, b, c) {
    const u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]], v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    const n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    const length = Math.hypot(...n) || 1;
    return [n[0] / length, n[1] / length, n[2] / length];
}

function flatTriangle(a, b, c) {
    const n = faceNormal(a, b, c);
    return {points: [a, b, c], normals: [n, n, n]};
}

/* A box centred on the origin. */
export function box(sx, sy, sz) {
    const x = sx / 2, y = sy / 2, z = sz / 2;
    const corners = [[-x, -y, -z], [x, -y, -z], [x, y, -z], [-x, y, -z], [-x, -y, z], [x, -y, z], [x, y, z], [-x, y, z]];
    const faces = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [2, 3, 7, 6], [1, 2, 6, 5], [3, 0, 4, 7]];
    const triangles = [];
    for (const [a, b, c, d] of faces) {
        triangles.push(flatTriangle(corners[a], corners[b], corners[c]), flatTriangle(corners[a], corners[c], corners[d]));
    }
    return build(triangles);
}

/*
A solid of revolution about the Z axis, from a profile of [radius, z]
points listed from bottom to top. Used for cylinders, cones and discs.
*/
export function lathe(profile, segments = 24) {
    const triangles = [];
    const point = (radius, z, k) => {
        const angle = 2 * Math.PI * k / segments;
        return [radius * Math.cos(angle), radius * Math.sin(angle), z];
    };
    for (let i = 0; i + 1 < profile.length; i++) {
        const [r0, z0] = profile[i], [r1, z1] = profile[i + 1];
        /* The outward normal of this band of the surface, in the (radius, z) plane. */
        const length = Math.hypot(r1 - r0, z1 - z0) || 1;
        const nr = (z1 - z0) / length, nz = -(r1 - r0) / length;
        for (let k = 0; k < segments; k++) {
            const normal = j => {
                const angle = 2 * Math.PI * j / segments;
                return [nr * Math.cos(angle), nr * Math.sin(angle), nz];
            };
            const a = point(r0, z0, k), b = point(r0, z0, k + 1), c = point(r1, z1, k + 1), d = point(r1, z1, k);
            if (r0 > 0) {
                triangles.push({points: [a, b, c], normals: [normal(k), normal(k + 1), normal(k + 1)]});
            }
            if (r1 > 0) {
                triangles.push({points: [a, c, d], normals: [normal(k), normal(k + 1), normal(k)]});
            }
        }
    }
    return build(triangles);
}

/* A closed cylinder along Z, from z = 0 to z = height. */
export function cylinder(radius, height, segments = 24) {
    return lathe([[0, 0], [radius, 0], [radius, height], [0, height]], segments);
}

/* A ball. */
export function sphere(radius, segments = 20) {
    const profile = [];
    for (let i = 0; i <= segments / 2; i++) {
        const angle = -Math.PI / 2 + Math.PI * i / (segments / 2);
        profile.push([radius * Math.cos(angle), radius * Math.sin(angle)]);
    }
    profile[0][0] = 0;
    profile[profile.length - 1][0] = 0;
    return lathe(profile, segments);
}

/* An arrow along +Z of unit length; scale it to the length wanted. */
export function arrowShape(shaft = 0.018, head = 0.05, headLength = 0.16, segments = 16) {
    return lathe([[0, 0], [shaft, 0], [shaft, 1 - headLength], [head, 1 - headLength], [0, 1]], segments);
}

/* Move, turn or stretch a shape by a 4x4 matrix (as in math.js). */
export function transform(shape, matrix) {
    const positions = new Float32Array(shape.positions.length);
    const normals = new Float32Array(shape.normals.length);
    for (let i = 0; i < shape.positions.length; i += 3) {
        const x = shape.positions[i], y = shape.positions[i + 1], z = shape.positions[i + 2];
        positions[i] = matrix[0] * x + matrix[4] * y + matrix[8] * z + matrix[12];
        positions[i + 1] = matrix[1] * x + matrix[5] * y + matrix[9] * z + matrix[13];
        positions[i + 2] = matrix[2] * x + matrix[6] * y + matrix[10] * z + matrix[14];
        const nx = shape.normals[i], ny = shape.normals[i + 1], nz = shape.normals[i + 2];
        const tx = matrix[0] * nx + matrix[4] * ny + matrix[8] * nz;
        const ty = matrix[1] * nx + matrix[5] * ny + matrix[9] * nz;
        const tz = matrix[2] * nx + matrix[6] * ny + matrix[10] * nz;
        const length = Math.hypot(tx, ty, tz) || 1;
        normals[i] = tx / length;
        normals[i + 1] = ty / length;
        normals[i + 2] = tz / length;
    }
    return {positions, normals};
}

/* Join shapes into one. */
export function merge(shapes) {
    const total = shapes.reduce((sum, shape) => sum + shape.positions.length, 0);
    const positions = new Float32Array(total), normals = new Float32Array(total);
    let offset = 0;
    for (const shape of shapes) {
        positions.set(shape.positions, offset);
        normals.set(shape.normals, offset);
        offset += shape.positions.length;
    }
    return {positions, normals};
}
