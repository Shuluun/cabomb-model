import * as THREE from 'three';

// All morphology is display-only: no measurements, drug positions, or ODE inputs.
function random(seed) {
  let state = seed >>> 0;
  return () => { state ^= state << 13; state ^= state >>> 17; state ^= state << 5; return (state >>> 0) / 4294967296; };
}

function noise(x, y, seed) {
  const hash = (a, b) => {
    let n = Math.imul(a + seed, 374761393) ^ Math.imul(b, 668265263);
    n = Math.imul(n ^ (n >>> 13), 1274126177);
    return ((n ^ (n >>> 16)) >>> 0) / 4294967296;
  };
  const ix = Math.floor(x), iy = Math.floor(y);
  let fx = x - ix, fy = y - iy;
  fx = fx * fx * (3 - 2 * fx); fy = fy * fy * (3 - 2 * fy);
  const a = hash(ix, iy) * (1 - fx) + hash(ix + 1, iy) * fx;
  const b = hash(ix, iy + 1) * (1 - fx) + hash(ix + 1, iy + 1) * fx;
  return a * (1 - fy) + b * fy;
}

function texture(kind, seed) {
  const size = 512;
  const canvas = document.createElement('canvas');
  canvas.width = canvas.height = size;
  const ctx = canvas.getContext('2d');
  const pixels = ctx.createImageData(size, size);
  const rng = random(seed);
  const palette = kind === 'tumor' ? [[198, 180, 164], [132, 87, 88]] : kind === 'nucleus' ? [[166, 144, 155], [105, 87, 106]] : [[181, 183, 170], [144, 150, 140]];
  for (let y = 0; y < size; y++) for (let x = 0; x < size; x++) {
    const u = x / size, v = y / size;
    const broad = noise(u * 6, v * 6, seed);
    const small = noise(u * 28 + broad * 2, v * 28, seed + 1);
    const fine = noise(u * 133, v * 133, seed + 2);
    let blend = kind === 'tumor' ? Math.max(0, Math.min(1, (broad - .43) * 2.5 + small * .22)) : broad * .35 + small * .4 + fine * .2;
    if (kind === 'nucleus') blend = Math.min(1, blend + Math.pow(Math.abs(small - .5) * 2, .65) * .3);
    const grain = (rng() - .5) * 17 + (fine - .5) * 19;
    const p = (y * size + x) * 4;
    for (let c = 0; c < 3; c++) pixels.data[p + c] = palette[0][c] * (1 - blend) + palette[1][c] * blend + grain;
    pixels.data[p + 3] = 255;
  }
  ctx.putImageData(pixels, 0, 0);
  // Fine, curved connective texture, rather than painted spherical lobules.
  const count = kind === 'tumor' ? 1250 : 500;
  for (let j = 0; j < count; j++) {
    const x = rng() * size, y = rng() * size, length = 5 + rng() * 85;
    const angle = kind === 'tumor' ? .8 * Math.sin(y / 60) + rng() * 1.8 : rng() * Math.PI * 2;
    ctx.beginPath(); ctx.moveTo(x, y);
    ctx.bezierCurveTo(x + length * .3 * Math.cos(angle + .3), y + length * .3 * Math.sin(angle + .3), x + length * .7 * Math.cos(angle - .2), y + length * .7 * Math.sin(angle - .2), x + length * Math.cos(angle), y + length * Math.sin(angle));
    ctx.strokeStyle = j % 4 ? 'rgba(235,220,200,.10)' : 'rgba(80,57,62,.12)';
    ctx.lineWidth = .35 + rng() * 1.2; ctx.stroke();
  }
  const result = new THREE.CanvasTexture(canvas);
  result.colorSpace = THREE.SRGBColorSpace;
  return result;
}

function material(color, extra = {}) {
  return new THREE.MeshStandardMaterial({ color, roughness: .82, metalness: 0, ...extra });
}

function radialOutline(angle, kind) {
  if (kind === 'tumor') return 1 + .095 * Math.sin(3 * angle + .5) + .055 * Math.cos(5 * angle - .7) + .025 * Math.sin(11 * angle);
  return 1 + .05 * Math.sin(3 * angle - .6) + .03 * Math.cos(5 * angle + .3) + .012 * Math.sin(9 * angle);
}

// A closed, irregular specimen with a front cut surface and a rounded back.
function specimen(kind, sx, sy, depth, faceZ) {
  const segments = 160, rings = 36;
  const positions = [], uvs = [], indices = [];
  const add = (x, y, z) => { positions.push(x, y, z); uvs.push(x / (sx * 2.5) + .5, y / (sy * 2.5) + .5); };
  for (let r = 0; r <= rings; r++) {
    const theta = r / rings * Math.PI * .62;
    const radius = Math.sin(theta);
    const z = -depth * Math.cos(theta);
    for (let s = 0; s <= segments; s++) {
      const a = s / segments * Math.PI * 2;
      const edge = radialOutline(a, kind);
      const ripple = .002 * Math.sin(a * 19 + r * .7) * Math.sin(theta);
      add(sx * radius * edge * Math.cos(a), sy * radius * edge * Math.sin(a), z + ripple);
    }
  }
  for (let r = 0; r < rings; r++) for (let s = 0; s < segments; s++) {
    const a = r * (segments + 1) + s, b = a + segments + 1;
    indices.push(a, a + 1, b, a + 1, b + 1, b);
  }
  const body = new THREE.BufferGeometry();
  body.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
  body.setAttribute('uv', new THREE.Float32BufferAttribute(uvs, 2)); body.setIndex(indices); body.computeVertexNormals();

  const facePositions = [], faceUvs = [], faceIndices = [];
  const rimRadius = Math.sin(Math.PI * .62);
  for (let r = 0; r <= rings; r++) for (let s = 0; s <= segments; s++) {
    const a = s / segments * Math.PI * 2, fraction = r / rings, radius = fraction * rimRadius;
    const edge = radialOutline(a, kind);
    const x = sx * radius * edge * Math.cos(a), y = sy * radius * edge * Math.sin(a);
    // Small uneven cut relief. Fixed in time; it does not encode amount.
    const ripple = .002 * Math.sin(a * 19 + rings * .7) * fraction ** 8 * rimRadius;
    const relief = (kind === 'tumor' ? .009 * (noise(x * 14, y * 14, 713) - .5) + .006 * Math.sin(x * 10 + y * 7) * Math.cos(y * 13) : .003 * Math.cos(x * 8 - y * 5)) * (1 - fraction ** 4);
    facePositions.push(x, y, faceZ + ripple + relief);
    faceUvs.push(x / (sx * 2.5) + .5, y / (sy * 2.5) + .5);
  }
  for (let r = 0; r < rings; r++) for (let s = 0; s < segments; s++) {
    const a = r * (segments + 1) + s, b = a + segments + 1;
    faceIndices.push(a, b, a + 1, a + 1, b, b + 1);
  }
  const face = new THREE.BufferGeometry();
  face.setAttribute('position', new THREE.Float32BufferAttribute(facePositions, 3));
  face.setAttribute('uv', new THREE.Float32BufferAttribute(faceUvs, 2)); face.setIndex(faceIndices); face.computeVertexNormals();
  return { body, face, rimRadius };
}

function tube(points, radius, mat, segments = 60) {
  return new THREE.Mesh(new THREE.TubeGeometry(new THREE.CatmullRomCurve3(points), segments, radius, 7, false), mat);
}

export function buildTumor() {
  const root = new THREE.Group();
  const tex = texture('tumor', 3717);
  const depth = .64, faceZ = -depth * Math.cos(Math.PI * .62);
  const geo = specimen('tumor', 1.16, .88, depth, faceZ);
  root.add(new THREE.Mesh(geo.body, material('#a68c83', { map: tex, bumpMap: tex, bumpScale: .06, side: THREE.DoubleSide })));
  root.add(new THREE.Mesh(geo.face, material('#eee6df', { map: tex, bumpMap: tex, bumpScale: .025, side: THREE.DoubleSide })));
  root.userData.signalMaterials = root.children.slice(0, 2).map(mesh => mesh.material);
  // Fine, irregular fibrous margin follows the tissue outline.
  const edge = [];
  for (let j = 0; j <= 180; j++) {
    const a = j / 180 * Math.PI * 2, r = geo.rimRadius * radialOutline(a, 'tumor');
    edge.push(new THREE.Vector3(1.16 * r * Math.cos(a), .88 * r * Math.sin(a), faceZ + .002 * Math.sin(a * 19 + 36 * .7) * geo.rimRadius));
  }
  root.add(tube(edge, .012, material('#b7a293'), 180));
  root.rotation.set(.24, -.36, -.14);
  return root;
}

function ribbon(points, widths, mat) {
  const verts = [], uv = [], indices = [];
  for (let i = 0; i < points.length; i++) {
    const tangent = points[Math.min(i + 1, points.length - 1)].clone().sub(points[Math.max(0, i - 1)]).normalize();
    const across = new THREE.Vector3(-tangent.y, tangent.x, .25 * Math.sin(i * .7)).normalize().multiplyScalar(widths[i] / 2);
    for (const side of [-1, 1]) { const p = points[i].clone().addScaledVector(across, side); verts.push(p.x, p.y, p.z); uv.push(i / (points.length - 1), (side + 1) / 2); }
    if (i) { const a = i * 2; indices.push(a - 2, a - 1, a, a - 1, a + 1, a); }
  }
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.Float32BufferAttribute(verts, 3)); geo.setAttribute('uv', new THREE.Float32BufferAttribute(uv, 2)); geo.setIndex(indices); geo.computeVertexNormals();
  return new THREE.Mesh(geo, mat);
}

function nucleusGeometry() {
  const geo = new THREE.SphereGeometry(1, 96, 64), p = geo.attributes.position;
  for (let i = 0; i < p.count; i++) {
    const x = p.getX(i), y = p.getY(i), z = p.getZ(i);
    const fold = 1 + .035 * Math.sin(x * 17 + y * 9) * Math.cos(y * 13 - z * 9) + .035 * Math.sin(y * 6 + x * 2);
    p.setXYZ(i, x * .40 * fold, y * .35 * fold, z * .23 * fold);
  }
  geo.computeVertexNormals(); return geo;
}

export function buildCell() {
  const root = new THREE.Group();
  const tex = texture('cytoplasm', 7113), nuclearTex = texture('nucleus', 8417);
  const depth = .64, faceZ = -depth * Math.cos(Math.PI * .62);
  const geo = specimen('cell', 1.12, .92, depth, faceZ);
  root.add(new THREE.Mesh(geo.body, material('#929c94', { map: tex, bumpMap: tex, bumpScale: .025, side: THREE.DoubleSide })));
  root.add(new THREE.Mesh(geo.face, material('#c8c7b8', { map: tex, bumpMap: tex, bumpScale: .012, side: THREE.DoubleSide })));
  root.userData.signalMaterials = root.children.slice(0, 2).map(mesh => mesh.material);
  const membrane = material('#8b9588', { roughness: .74 });
  for (const inset of [1, .984]) {
    const edge = [];
    for (let j = 0; j <= 180; j++) {
      const a = j / 180 * Math.PI * 2, r = inset * geo.rimRadius * radialOutline(a, 'cell');
      edge.push(new THREE.Vector3(1.12 * r * Math.cos(a), .92 * r * Math.sin(a), faceZ + .002 * Math.sin(a * 19 + 36 * .7) * geo.rimRadius));
    }
    root.add(tube(edge, inset === 1 ? .017 : .006, membrane, 180));
  }
  const nucleus = new THREE.Mesh(nucleusGeometry(), material('#d2c2ce', { map: nuclearTex, bumpMap: nuclearTex, bumpScale: .023 }));
  nucleus.position.set(-.22, .025, faceZ + .11); nucleus.rotation.z = -.16; root.add(nucleus);
  const nucleolus = new THREE.Mesh(new THREE.SphereGeometry(.088, 32, 24), material('#675767', { map: nuclearTex, bumpMap: nuclearTex, bumpScale: .01 }));
  nucleolus.scale.set(1.1, .85, .6); nucleolus.position.set(-.27, .065, faceZ + .326); root.add(nucleolus);

  const erMat = material('#bbb2aa', { map: tex, bumpMap: tex, bumpScale: .014, side: THREE.DoubleSide });
  const erEdge = material('#8a8d85', { roughness: .9 });
  for (let j = 0; j < 8; j++) {
    const points = [], widths = [];
    for (let k = 0; k <= 64; k++) {
      const a = -.8 + k / 64 * 4.95;
      const r = .44 + j * .043 + .022 * Math.sin(a * 7 + j * .7);
      points.push(new THREE.Vector3(-.22 + r * Math.cos(a), .025 + r * .73 * Math.sin(a), faceZ + .028 + .025 * Math.sin(a * 5 + j * .5)));
      widths.push(.037 + .012 * Math.cos(a * 6 + j));
    }
    root.add(ribbon(points, widths, erMat)); root.add(tube(points, .0035, erEdge, 64));
  }
  const rng = random(9187), mitoMat = material('#9e8a6e', { map: tex, bumpMap: tex, bumpScale: .012 }), cristaeMat = material('#655f53');
  const mitochondria = [[.61,.13,-.3],[.59,-.28,.75],[.27,-.54,1.4],[-.3,-.57,-1.1],[-.81,-.17,.4],[-.76,.28,.8],[.1,.53,-.7],[.51,.42,1.2],[.87,-.05,.2],[-.43,.51,.2]];
  for (const [x, y, angle] of mitochondria) {
    const group = new THREE.Group(), path = [];
    const length = .12 + rng() * .08;
    for (let k = 0; k <= 20; k++) { const t = k / 20; path.push(new THREE.Vector3((t - .5) * length, .022 * Math.sin(t * Math.PI * 1.5), 0)); }
    group.add(tube(path, .038, mitoMat, 24));
    for (let k = 0; k < 6; k++) {
      const t = (k + .5) / 6, center = path[Math.round(t * 20)], fold = [];
      for (let j = 0; j < 7; j++) fold.push(new THREE.Vector3(center.x + .006 * Math.sin(j / 6 * Math.PI * 2), center.y + (j / 6 - .5) * .058, .032 + .004 * Math.sin(j)));
      group.add(tube(fold, .003, cristaeMat, 10));
    }
    group.position.set(x, y, faceZ + .05); group.rotation.z = angle; root.add(group);
  }
  const golgiMat = material('#b4a994', { side: THREE.DoubleSide, bumpMap: tex, bumpScale: .008 });
  for (let j = 0; j < 5; j++) {
    const curve = [], widths = [];
    for (let k = 0; k <= 30; k++) { const t = k / 30; curve.push(new THREE.Vector3(.29 + t * .33 + .012 * Math.sin(t * 14 + j), -.11 - j * .035 + .055 * Math.sin(t * Math.PI), faceZ + .04 + .02 * Math.sin(t * Math.PI))); widths.push(.017 + .008 * Math.sin(t * Math.PI)); }
    root.add(ribbon(curve, widths, golgiMat));
  }
  root.rotation.set(.25, -.38, -.1);
  return root;
}
