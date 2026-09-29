// 3-D view: transparent isosurface of the binned volume plus the three current
// slices as textured planes and the line cut, all in the lattice geometry. Geometry is built in
// display coordinates (r.l.u. along each axis) inside a group whose matrix maps
// them to Cartesian reciprocal space, so the isosurface and slices are sheared
// exactly like the 2-D plots. Rendering happens on demand.

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

// Outline colors follow the page's axis hues: H amber, K blue, L green.
const PLANE_COLORS = [0xd98a0b, 0x2f74e6, 0x13a36b];
// The line cut: a color apart from the axes, the isosurface and the colormaps.
const CUT_COLOR = 0xc2185b;

export class View3D {
  constructor(canvas) {
    this.canvas = canvas;
    this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true, preserveDrawingBuffer: true });
    this.renderer.setPixelRatio(devicePixelRatio || 1);
    this.renderer.setClearColor(0xffffff, 1);
    this.renderer.localClippingEnabled = true;
    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(35, 1, 0.01, 100);
    this.camera.up.set(0, 0, 1);
    this.scene.add(this.camera);
    this.scene.add(new THREE.AmbientLight(0xffffff, 1.2));
    const light = new THREE.DirectionalLight(0xffffff, 2.2);
    light.position.set(1, 1.5, 2);
    this.camera.add(light);
    this.world = new THREE.Group();
    this.world.matrixAutoUpdate = false;
    this.scene.add(this.world);
    this.labels = new THREE.Group();
    this.scene.add(this.labels);
    this.controls = new OrbitControls(this.camera, canvas);
    this.controls.addEventListener('change', () => this.render());
    this.slices = new Map();
    this.sliceOpacity = 1;
    this.mesh = null;
    this.frame = null;
    this.cut = null;
    new ResizeObserver(() => this.resize()).observe(canvas);
  }

  /**
   * Display->Cartesian matrix T (row-major 3x3), the view box in display
   * coordinates [[lo, hi] x 3] and axis labels. Re-frames the camera when the
   * box or basis changes.
   */
  setFrame(T, box, labels) {
    const key = JSON.stringify([T, box, labels]);
    if (key === this.frameKey) return;
    // Re-frame the camera only when the box changes, not for small basis changes.
    const reframe = JSON.stringify(box) !== this.boxKey;
    this.frameKey = key;
    this.boxKey = JSON.stringify(box);
    this.T = T;
    this.box = box;
    this.world.matrix.set(T[0], T[1], T[2], 0, T[3], T[4], T[5], 0, T[6], T[7], T[8], 0, 0, 0, 0, 1);
    this.world.matrixWorldNeedsUpdate = true;

    if (this.frame) { this.world.remove(this.frame); this.frame.geometry.dispose(); }
    const corner = (i) => box.map((r, d) => r[(i >> d) & 1]);
    const points = [];
    for (let i = 0; i < 8; i++) {
      for (const bit of [1, 2, 4]) if (!(i & bit)) points.push(...corner(i), ...corner(i | bit));
    }
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(points, 3));
    this.frame = new THREE.LineSegments(geometry, new THREE.LineBasicMaterial({ color: 0x9aa6b2 }));
    this.world.add(this.frame);

    // Clip the isosurface to the view box: x_d = g_d . P with g_d the rows of T^-1.
    const inv = new THREE.Matrix3().set(...T).invert().elements; // column-major
    this.clipping = box.flatMap(([lo, hi], d) => {
      const g = new THREE.Vector3(inv[d], inv[d + 3], inv[d + 6]), len = g.length();
      return [new THREE.Plane(g.clone().divideScalar(len), -lo / len), new THREE.Plane(g.clone().divideScalar(-len), hi / len)];
    });
    if (this.mesh) this.mesh.material.clippingPlanes = this.clipping;
    for (const part of this.cutParts()) part.material.clippingPlanes = this.clipping;

    const center = this.toCartesian(box.map(([lo, hi]) => (lo + hi) / 2));
    this.radius = Math.max(...Array.from({ length: 8 }, (_, i) => this.toCartesian(corner(i)).distanceTo(center)));
    for (const sprite of [...this.labels.children]) { sprite.material.map.dispose(); sprite.material.dispose(); }
    this.labels.clear();
    // Label each axis at the middle of its box edge through the low corner,
    // pushed outward from the box center (like matplotlib's 3-D axes).
    labels.forEach((text, d) => {
      const at = box.map((r) => r[0]);
      at[d] = (box[d][0] + box[d][1]) / 2;
      const edge = this.toCartesian(at);
      const sprite = textSprite(text, this.radius * 0.045, PLANE_COLORS[d]);
      sprite.position.copy(edge).add(edge.clone().sub(center).normalize().multiplyScalar(this.radius * 0.1));
      this.labels.add(sprite);
    });
    if (reframe) this.resetView();
    else this.render();
  }

  toCartesian(p) {
    const T = this.T;
    return new THREE.Vector3(
      T[0] * p[0] + T[1] * p[1] + T[2] * p[2],
      T[3] * p[0] + T[4] * p[1] + T[5] * p[2],
      T[6] * p[0] + T[7] * p[1] + T[8] * p[2]);
  }

  resetView() {
    if (!this.box) return;
    const center = this.toCartesian(this.box.map(([lo, hi]) => (lo + hi) / 2));
    const distance = this.radius / Math.sin((this.camera.fov / 2) * Math.PI / 180) * 0.8;
    this.controls.target.copy(center);
    this.camera.position.copy(center).add(new THREE.Vector3(1.1, -1.6, 0.9).normalize().multiplyScalar(distance));
    this.camera.near = distance / 100;
    this.camera.far = distance * 10;
    this.camera.updateProjectionMatrix();
    this.controls.update();
    this.render();
  }

  /** Isosurface in display coordinates; null clears it. */
  setMesh(positions, indices, opacity) {
    if (this.mesh) { this.world.remove(this.mesh); this.mesh.geometry.dispose(); this.mesh.material.dispose(); this.mesh = null; }
    if (positions?.length) {
      const geometry = new THREE.BufferGeometry();
      geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
      geometry.setIndex(new THREE.BufferAttribute(indices, 1));
      geometry.computeVertexNormals();
      const material = new THREE.MeshPhongMaterial({
        color: 0x2f6fd0, specular: 0x444444, shininess: 50, side: THREE.DoubleSide,
        transparent: true, opacity, depthWrite: false, clippingPlanes: this.clipping ?? [],
      });
      this.mesh = new THREE.Mesh(geometry, material);
      this.world.add(this.mesh);
    }
    this.render();
  }

  /**
   * Slice planes are opaque at 1; below that they blend and stop writing depth,
   * so the isosurface shows through. No-data pixels (texture alpha 0) stay cut
   * out because the alpha test sits below the plane opacity.
   */
  setSliceOpacity(opacity) {
    this.sliceOpacity = opacity;
    for (const { plane } of this.slices.values()) styleSlice(plane.material, opacity);
    this.render();
  }

  setOpacity(opacity) {
    if (this.mesh) this.mesh.material.opacity = opacity;
    this.render();
  }

  /**
   * Slice planes: {fixed, x, y, center, image, key, u: [lo, hi], v: [lo, hi],
   * ex: [first, last edge], ey: [...]}. Pixels with zero alpha (no data) are cut out.
   */
  setSlices(list, visible) {
    const keep = new Set();
    for (const s of list) {
      keep.add(s.fixed);
      let entry = this.slices.get(s.fixed);
      if (!entry) {
        const material = new THREE.MeshBasicMaterial({ side: THREE.DoubleSide });
        styleSlice(material, this.sliceOpacity);
        const plane = new THREE.Mesh(new THREE.BufferGeometry(), material);
        const outline = new THREE.LineLoop(new THREE.BufferGeometry(), new THREE.LineBasicMaterial({ color: PLANE_COLORS[s.fixed] }));
        entry = { plane, outline };
        this.world.add(plane, outline);
        this.slices.set(s.fixed, entry);
      }
      if (entry.key !== s.key) {
        entry.key = s.key;
        entry.plane.material.map?.dispose();
        const texture = new THREE.CanvasTexture(s.image);
        texture.flipY = false;
        texture.magFilter = THREE.NearestFilter;
        texture.minFilter = THREE.LinearFilter;
        texture.generateMipmaps = false;
        texture.colorSpace = THREE.SRGBColorSpace;
        entry.plane.material.map = texture;
        entry.plane.material.needsUpdate = true;
      }
      const at = (u, v) => { const p = [0, 0, 0]; p[s.x] = u; p[s.y] = v; p[s.fixed] = s.center; return p; };
      const corners = [[s.u[0], s.v[0]], [s.u[1], s.v[0]], [s.u[1], s.v[1]], [s.u[0], s.v[1]]];
      const positions = corners.flatMap(([u, v]) => at(u, v));
      const uv = corners.flatMap(([u, v]) => [(u - s.ex[0]) / (s.ex[1] - s.ex[0]), (v - s.ey[0]) / (s.ey[1] - s.ey[0])]);
      entry.plane.geometry.dispose();
      entry.plane.geometry = new THREE.BufferGeometry();
      entry.plane.geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
      entry.plane.geometry.setAttribute('uv', new THREE.Float32BufferAttribute(uv, 2));
      entry.plane.geometry.setIndex([0, 1, 2, 0, 2, 3]);
      entry.outline.geometry.dispose();
      entry.outline.geometry = new THREE.BufferGeometry();
      entry.outline.geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
      entry.plane.visible = entry.outline.visible = visible;
    }
    for (const [fixed, entry] of this.slices) {
      if (!keep.has(fixed)) entry.plane.visible = entry.outline.visible = false;
    }
    this.render();
  }

  /**
   * The line cut: `line`, its two ends, `rod`, the ends of the rod of voxels it
   * averages (its first and last bin edges on the line), and `pierce`, where it
   * passes through the other slices ({ at, axis }: a point, and the display axis
   * normal to that slice), all in display coordinates, and `radius`, the rod's
   * radius in Cartesian units; null hides it.
   */
  setCut(cut) {
    if (cut && !this.T) return; // no frame yet
    if (!this.cut && cut) {
      const shape = (geometry, material, order = 0) => {
        material.clippingPlanes = this.clipping ?? [];
        const mesh = new THREE.Mesh(geometry, material);
        mesh.renderOrder = order;
        return mesh;
      };
      // The rod and its ends are drawn twice, so they read as passing through the
      // slices: shaded where nothing is in front of them, and as a faint ghost,
      // drawn last without the depth test, where a slice hides them.
      const solid = () => new THREE.MeshPhongMaterial({ color: CUT_COLOR, specular: 0x444444, shininess: 40 });
      const ghost = () => new THREE.MeshBasicMaterial({ color: CUT_COLOR, transparent: true, opacity: 0.28, depthTest: false, depthWrite: false });
      // Built in Cartesian space, so they stay round under the lattice shear.
      const rod = new THREE.CylinderGeometry(1, 1, 1, 20), end = new THREE.SphereGeometry(1, 20, 14);
      // The rod of voxels it averages: an open translucent tube with its end circles outlined.
      const tube = new THREE.CylinderGeometry(1, 1, 1, 48, 1, true);
      this.cut = {
        sleeve: shape(tube, new THREE.MeshBasicMaterial({ color: CUT_COLOR, transparent: true, opacity: 0.22, depthWrite: false, side: THREE.DoubleSide })),
        rims: new THREE.LineSegments(new THREE.EdgesGeometry(tube, 15), new THREE.LineBasicMaterial({ color: CUT_COLOR, transparent: true, opacity: 0.6, clippingPlanes: this.clipping ?? [] })),
        rods: [shape(rod, solid()), shape(rod, ghost(), 10)],
        ends: [shape(end, solid()), shape(end, solid()), shape(end, ghost(), 10), shape(end, ghost(), 10)],
        // A collar where the rod passes through another slice, lifted off the plane so it does not flicker.
        collars: [0, 1].map(() => shape(new THREE.RingGeometry(1.05, 2.3, 32), new THREE.MeshBasicMaterial({
          color: 0xffffff, side: THREE.DoubleSide, polygonOffset: true, polygonOffsetFactor: -4, polygonOffsetUnits: -4,
        }))),
      };
      this.scene.add(this.cut.sleeve, this.cut.rims, ...this.cut.rods, ...this.cut.ends, ...this.cut.collars);
    }
    for (const part of this.cutParts()) part.visible = !!cut;
    if (cut) {
      const { sleeve, rims, rods, ends } = this.cut;
      // A unit cylinder along y, stretched between two Cartesian points with the given radius.
      const span = (mesh, from, to, radius) => {
        mesh.position.copy(from).add(to).multiplyScalar(0.5);
        mesh.scale.set(radius, from.distanceTo(to), radius);
        mesh.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), to.clone().sub(from).normalize());
      };
      const [a, b] = cut.line.map((p) => this.toCartesian(p)), [s0, s1] = cut.rod.map((p) => this.toCartesian(p));
      const r = this.radius * 0.006;
      for (const part of [sleeve, rims]) span(part, s0, s1, cut.radius);
      for (const rod of rods) span(rod, a, b, r);
      ends.forEach((end, i) => { end.position.copy(i % 2 ? b : a); end.scale.setScalar(2.2 * r); });
      // Each collar lies in its slice: its normal is the gradient of that display coordinate.
      const inv = new THREE.Matrix3().set(...this.T).invert().elements; // column-major
      this.cut.collars.forEach((collar, i) => {
        const through = cut.pierce[i];
        collar.visible = !!through;
        if (!through) return;
        const d = through.axis, normal = new THREE.Vector3(inv[d], inv[d + 3], inv[d + 6]).normalize();
        collar.position.copy(this.toCartesian(through.at));
        collar.quaternion.setFromUnitVectors(new THREE.Vector3(0, 0, 1), normal);
        collar.scale.setScalar(r);
      });
    }
    this.render();
  }

  cutParts() {
    return this.cut ? [this.cut.sleeve, this.cut.rims, ...this.cut.rods, ...this.cut.ends, ...this.cut.collars] : [];
  }

  resize() {
    const w = this.canvas.clientWidth, h = this.canvas.clientHeight;
    if (!w || !h) return;
    this.renderer.setSize(w, h, false);
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
    this.render();
  }

  render() {
    if (this.pending) return;
    this.pending = requestAnimationFrame(() => {
      this.pending = null;
      this.renderer.render(this.scene, this.camera);
    });
  }

  snapshot(callback) {
    this.renderer.render(this.scene, this.camera);
    this.canvas.toBlob(callback);
  }
}

function styleSlice(material, opacity) {
  const blend = opacity < 1;
  material.transparent = blend;
  material.opacity = opacity;
  material.depthWrite = !blend;
  material.alphaTest = blend ? opacity / 2 : 0.5;
  material.needsUpdate = true;
}

function textSprite(text, height, color) {
  const canvas = document.createElement('canvas');
  const ctx = canvas.getContext('2d');
  ctx.font = '600 44px system-ui, sans-serif';
  canvas.width = Math.ceil(ctx.measureText(text).width) + 16;
  canvas.height = 60;
  ctx.font = '600 44px system-ui, sans-serif';
  ctx.fillStyle = `#${color.toString(16).padStart(6, '0')}`;
  ctx.textBaseline = 'middle';
  ctx.fillText(text, 8, 32);
  const texture = new THREE.CanvasTexture(canvas);
  texture.colorSpace = THREE.SRGBColorSpace;
  const sprite = new THREE.Sprite(new THREE.SpriteMaterial({ map: texture, depthTest: false }));
  sprite.scale.set(height * canvas.width / canvas.height, height, 1);
  return sprite;
}
