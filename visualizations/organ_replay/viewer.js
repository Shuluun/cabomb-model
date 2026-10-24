import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { GLTFLoader } from 'three/addons/loaders/GLTFLoader.js';
import { buildTumor, buildCell } from './morphology.js';
import { attachAmountSignal, MAX_SIGNAL_BLEND } from './amount_signal.js';

const data = window.CABOMB_DATA;
const points = data.points;
const colors = ['#68d4df', '#68d4df', '#68d4df'];
const views = [];
const charts = [];
let currentTime = 24;
let playing = false;
let previousStamp = 0;
let dirty = true;
const svgNS = 'http://www.w3.org/2000/svg';
const playButton = document.getElementById('play');
const slider = document.getElementById('timeline');
const amountCeilings = data.displayCeilings;
const endTime = points[points.length - 1][0];
// Fixed display texture only; no amount or spatial uptake is encoded.
function displayRandom(seed) {
  let value = seed >>> 0;
  return () => { value ^= value << 13; value ^= value >>> 17; value ^= value << 5; return (value >>> 0) / 4294967296; };
}

function makeTissueTexture() {
  const canvas = document.createElement('canvas');
  canvas.width = canvas.height = 128;
  const ctx = canvas.getContext('2d');
  const image = ctx.createImageData(128, 128);
  const random = displayRandom(41987);
  for (let i = 0; i < image.data.length; i += 4) {
    const value = 125 + Math.round(random() * 40);
    image.data[i] = image.data[i + 1] = image.data[i + 2] = value;
    image.data[i + 3] = 255;
  }
  ctx.putImageData(image, 0, 0);
  const texture = new THREE.CanvasTexture(canvas);
  texture.wrapS = texture.wrapT = THREE.RepeatWrapping;
  texture.repeat.set(4, 4);
  return texture;
}

const tissueTexture = makeTissueTexture();
function sample(time) {
  let lo = 0, hi = points.length - 1;
  while (lo + 1 < hi) { const mid = (lo + hi) >> 1; if (points[mid][0] <= time) lo = mid; else hi = mid; }
  const a = points[lo], b = points[hi];
  const f = Math.max(0, Math.min(1, (time - a[0]) / (b[0] - a[0])));
  return a.slice(1).map((value, i) => value + f * (b[i + 1] - value));
}

function svgEl(tag, attrs, parent) {
  const el = document.createElementNS(svgNS, tag);
  for (const [key, value] of Object.entries(attrs)) el.setAttribute(key, value);
  if (parent) parent.append(el);
  return el;
}

function makeChart(index) {
  const svg = document.getElementById(`chart-${index}`);
  svg.setAttribute('viewBox', '0 0 360 130');
  const left = 36, right = 346, top = 10, bottom = 104;
  const ceiling = amountCeilings[index];
  const x = time => left + time / endTime * (right - left);
  const y = value => bottom - value / ceiling * (bottom - top);
  for (const fraction of [0, .5, 1]) {
    const yy = y(ceiling * fraction);
    svgEl('line', { x1: left, x2: right, y1: yy, y2: yy, class: 'grid' }, svg);
    const label = svgEl('text', { x: left - 8, y: yy + 3, 'text-anchor': 'end' }, svg);
    label.textContent = (ceiling * fraction).toLocaleString('en-US', { maximumFractionDigits: 2 });
  }
  for (const time of data.timeTicks) {
    const label = svgEl('text', { x: x(time), y: bottom + 18, 'text-anchor': 'middle' }, svg);
    label.textContent = time;
  }
  const curve = points.map((row, i) => `${i ? 'L' : 'M'}${x(row[0]).toFixed(2)},${y(row[index + 1]).toFixed(2)}`).join(' ');
  const defs = svgEl('defs', {}, svg);
  const gradient = svgEl('linearGradient', { id: `fill-${index}`, x1: 0, y1: 0, x2: 0, y2: 1 }, defs);
  svgEl('stop', { offset: '0%', 'stop-color': colors[index], 'stop-opacity': .29 }, gradient);
  svgEl('stop', { offset: '100%', 'stop-color': colors[index], 'stop-opacity': .015 }, gradient);
  const history = svgEl('rect', { x: left, y: 0, width: 0, height: 130 }, svgEl('clipPath', { id: `past-${index}` }, defs));
  svgEl('path', { d: curve, stroke: colors[index], class: 'curve future' }, svg);
  const past = svgEl('g', { 'clip-path': `url(#past-${index})` }, svg);
  svgEl('path', { d: `${curve} L${right},${bottom} L${left},${bottom} Z`, fill: `url(#fill-${index})` }, past);
  svgEl('path', { d: curve, stroke: colors[index], class: 'curve' }, past);
  const cursor = svgEl('line', { x1: x(currentTime), x2: x(currentTime), y1: top, y2: bottom, class: 'cursor' }, svg);
  const dot = svgEl('circle', { r: 4.5, fill: colors[index], class: 'dot' }, svg);
  charts.push({ cursor, dot, x, y, history });
  svg.addEventListener('pointerdown', event => {
    const rect = svg.getBoundingClientRect();
    setPlaying(false);
    const scale = Math.min(rect.width / 360, rect.height / 130);
    const inset = (rect.width - 360 * scale) / 2;
    setTime(Math.max(0, Math.min(endTime, (((event.clientX - rect.left - inset) / scale) - left) / (right - left) * endTime)));
  });
}

function makeView(index) {
  const container = document.getElementById(`view-${index}`);
  let renderer;
  try {
    renderer = new THREE.WebGLRenderer({ alpha: true, antialias: true, preserveDrawingBuffer: true });
  } catch (error) {
    const fallback = document.createElement('div');
    fallback.className = 'fallback';
    fallback.textContent = '当前设备未能显示三维模型；下方数值与曲线仍可查看。';
    container.append(fallback);
    console.warn('WebGL view unavailable', index, error.message);
    return null;
  }
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 1.75));
  renderer.setClearColor(0x000000, 0);
  renderer.outputColorSpace = THREE.SRGBColorSpace;
  renderer.toneMapping = THREE.ACESFilmicToneMapping;
  renderer.toneMappingExposure = 1.05;
  renderer.domElement.setAttribute('aria-hidden', 'true');
  container.append(renderer.domElement);
  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(32, 1, .1, 50);
  camera.position.set(0, .2, 3.8);
  scene.add(new THREE.HemisphereLight(0xcde8ff, 0x4d3834, 1.25));
  const key = new THREE.DirectionalLight(0xffeddf, 2.1);
  key.position.set(-3, 5, 5);
  scene.add(key);
  const rim = new THREE.DirectionalLight(0xb8d0d8, 1.6);
  rim.position.set(4, 2, -3);
  scene.add(rim);
  const controls = new OrbitControls(camera, renderer.domElement);
  controls.enablePan = false;
  controls.enableZoom = false;
  controls.minDistance = 3.1;
  controls.maxDistance = 7;
  controls.target.set(0, 0, 0);
  controls.enableDamping = false;
  controls.addEventListener('change', () => { dirty = true; });
  const view = { renderer, scene, camera, controls, root: null, signalBindings: [] };
  const resize = () => {
    const width = container.clientWidth, height = container.clientHeight;
    renderer.setSize(width, height, false);
    camera.aspect = width / height;
    camera.updateProjectionMatrix();
    dirty = true;
  };
  new ResizeObserver(resize).observe(container);
  resize();
  return view;
}

function attachRoot(view, root, index) {
  if (!view) return;
  view.root = root;
  const box = new THREE.Box3().setFromObject(root);
  const center = box.getCenter(new THREE.Vector3());
  const size = box.getSize(new THREE.Vector3());
  const scale = 2.15 / Math.max(size.x, size.y, size.z);
  root.position.sub(center);
  const wrapper = new THREE.Group();
  wrapper.add(root);
  wrapper.scale.setScalar(scale);
  wrapper.position.x = 0;
  view.scene.add(wrapper);
  view.signalBindings = attachAmountSignal(root);
  const container = document.getElementById(`view-${index}`);
  container.dataset.originalGeometry = 'true';
  container.dataset.signalMaterials = view.signalBindings.length;
  container.dataset.baseColor = view.signalBindings[0]?.material.color.getHexString() || '';
  container.dataset.surfaceOpacity = view.signalBindings[0]?.material.opacity ?? 1;
  const floorCanvas = document.createElement('canvas');
  floorCanvas.width = floorCanvas.height = 128;
  const ctx = floorCanvas.getContext('2d');
  const gradient = ctx.createRadialGradient(64, 64, 3, 64, 64, 64);
  gradient.addColorStop(0, 'rgba(0,0,0,.55)');
  gradient.addColorStop(1, 'rgba(0,0,0,0)');
  ctx.fillStyle = gradient; ctx.fillRect(0, 0, 128, 128);
  const shadow = new THREE.Mesh(new THREE.PlaneGeometry(3.2, 1.7), new THREE.MeshBasicMaterial({ map: new THREE.CanvasTexture(floorCanvas), transparent: true, depthWrite: false }));
  shadow.rotation.x = -Math.PI / 2;
  shadow.position.y = -1.13;
  view.scene.add(shadow);
  containerReady(index);
  setTime(currentTime);
}

let loaded = 0;
function containerReady(index) {
  document.getElementById(`view-${index}`).dataset.ready = 'true';
  document.getElementById('load-progress').style.width = `${++loaded / 3 * 100}%`;
  if (loaded === 3) { document.body.dataset.modelsReady = 'true'; document.getElementById('load-progress').style.opacity = 0; }
}

function setTime(time) {
  currentTime = Math.max(0, Math.min(endTime, Number(time)));
  slider.value = currentTime;
  slider.style.setProperty('--progress', `${currentTime / endTime * 100}%`);
  document.querySelectorAll('[data-time]').forEach(button => button.classList.toggle('active', Math.abs(Number(button.dataset.time) - currentTime) < .05));
  slider.setAttribute('aria-valuetext', `${currentTime.toFixed(1)}小时`);
  document.getElementById('time-value').textContent = currentTime.toFixed(1);
  const values = sample(currentTime);
  values.forEach((value, i) => {
    const formatted = value === 0 ? '0.000' : value.toFixed(3);
    document.getElementById(`amount-${i}`).textContent = formatted;
    document.getElementById(`amount-${i}`).dataset.exactValue = value;
    const c = charts[i];
    c.cursor.setAttribute('x1', c.x(currentTime)); c.cursor.setAttribute('x2', c.x(currentTime));
    c.dot.setAttribute('cx', c.x(currentTime)); c.dot.setAttribute('cy', c.y(value));
    document.getElementById(`chart-${i}`).setAttribute('aria-label', `${['肝关联载体', '肿瘤组织载体', '细胞内活性药物'][i]}，${currentTime.toFixed(1)}小时为${formatted}%`);
    const ratio = Math.max(0, Math.min(1, value / amountCeilings[i]));
    const container = document.getElementById(`view-${i}`);
    container.dataset.amountRatio = ratio;
    container.setAttribute('aria-label', `${['肝脏参考模型', '肿瘤组织切面示意', '细胞质剖面示意'][i]}，青蓝色叠色表示当前总体含量${formatted}%，量程0至${amountCeilings[i]}%`);
    const view = views[i];
    if (view) {
      view.signalBindings.forEach(({ amount }) => { amount.value = ratio; });
      container.dataset.signalRatio = ratio;
      container.dataset.maxSignalBlend = ratio * MAX_SIGNAL_BLEND;
      container.dataset.baseColor = view.signalBindings[0]?.material.color.getHexString() || '';
      container.dataset.surfaceOpacity = view.signalBindings[0]?.material.opacity ?? 1;
      container.dataset.surfacesOpaque = view.signalBindings.every(({ material }) => !material.transparent && material.opacity === 1);

    }
    c.history.setAttribute('width', Math.max(0, c.x(currentTime) - 36));
  });
  dirty = true;
}

function setPlaying(value) {
  playing = value;
  previousStamp = 0;
  playButton.querySelector('span').textContent = value ? '暂停' : '播放';
  playButton.setAttribute('aria-label', value ? '暂停时间回放' : '播放时间回放');
  playButton.querySelector('path').setAttribute('d', value ? 'M3 2h4v12H3zM10 2h4v12h-4z' : 'M4 2l10 6-10 6z');
}

for (let i = 0; i < 3; i++) { makeChart(i); views.push(makeView(i)); }
attachRoot(views[1], buildTumor(), 1);
attachRoot(views[2], buildCell(), 2);
try {
  const bytes = Uint8Array.from(atob(window.CABOMB_LIVER), char => char.charCodeAt(0));
  new GLTFLoader().parse(bytes.buffer, '', gltf => {
    const root = gltf.scene;
    const signalMaterials = [];
    root.traverse(object => {
      if (!object.isMesh) return;
      object.material = new THREE.MeshStandardMaterial({ color: '#96564a', roughness: .44, metalness: 0, bumpMap: tissueTexture, bumpScale: .025 });
      signalMaterials.push(object.material);
    });
    root.userData.signalMaterials = signalMaterials;
    root.rotation.set(0, .12, -.08);
    attachRoot(views[0], root, 0);
  }, error => { document.getElementById('view-0').textContent = '肝脏参考模型加载失败'; console.error(error); });
} catch (error) { console.error(error); }
setTime(currentTime);
slider.addEventListener('input', event => { setPlaying(false); setTime(event.target.value); });
document.querySelectorAll('[data-time]').forEach(button => button.addEventListener('click', () => { setPlaying(false); setTime(button.dataset.time); }));
playButton.addEventListener('click', () => { if (!playing && currentTime >= endTime) setTime(0); setPlaying(!playing); });
document.addEventListener('visibilitychange', () => { if (document.hidden) setPlaying(false); });
function frame(stamp) {
  if (playing) {
    if (previousStamp) setTime(currentTime + Math.min((stamp - previousStamp) / 1000, .1) * endTime / 30);
    previousStamp = stamp;
    if (currentTime >= endTime) setPlaying(false);
  }
  if (dirty) { views.forEach(view => { if (view) view.renderer.render(view.scene, view.camera); }); dirty = false; }
  requestAnimationFrame(frame);
}
requestAnimationFrame(frame);
