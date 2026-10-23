"""
Anopheles 4th-instar larva — parametric model for Blender (bpy, Blender 4.x/5.x)

Units: mm, real size (total length ~5.6 mm incl. anal papillae).
Axes: +X = posterior (head at -X), +Y = left/right, +Z = dorsal.
The larva rests parallel to the water surface (XY plane), dorsal side up.

Run:
  blender -b -P anopheles_larva.py -- --out ./out              (model + renders)
  blender -b -P anopheles_larva.py -- --out ./out --print       (thickened appendages, for printing)
or from Blender's Text Editor: set PARAMS below and Run Script.

Hairs (setae) are omitted except optional simplified palmate hairs and lateral setae.
"""
import bpy, bmesh, math, sys, os, argparse
from mathutils import Vector, Euler

# ---------------- parameters ----------------
PARAMS = dict(
    ring_verts=48,          # vertices per cross-section
    flatten_exp=2.4,        # superellipse exponent (2 = ellipse, larger = boxier)
    seg_len=0.46,           # abdominal segment I–VII length (mm)
    constriction=0.12,      # depth of intersegmental grooves (0–0.3)
    appendage_thick=1.0,    # radius multiplier for antennae, brushes, papillae (print: 2–3)
    palmate=True,           # simplified palmate hairs on III–VII (dorsal discs)
    lateral_setae=False,    # thin lateral hairs on thorax and abdomen
    print_scale=10.0,       # scale for the print STL
    voxel=0.010,            # voxel size (mm, real scale) for watertight union
)

# ---------------- helpers ----------------
def mat(name, rgb, rough=0.45, trans=0.0, alpha=1.0):
    m = bpy.data.materials.get(name) or bpy.data.materials.new(name)
    m.use_nodes = True
    b = m.node_tree.nodes.get("Principled BSDF")
    b.inputs["Base Color"].default_value = (*rgb, 1)
    b.inputs["Roughness"].default_value = rough
    for k in ("Transmission Weight", "Transmission"):
        if k in b.inputs:
            b.inputs[k].default_value = trans
            break
    b.inputs["Alpha"].default_value = alpha
    return m

def ellipsoid(name, center, radii, rot=(0, 0, 0), material=None, seg=24):
    bpy.ops.mesh.primitive_uv_sphere_add(segments=seg, ring_count=seg // 2, radius=1,
                                         location=center, rotation=rot)
    o = bpy.context.active_object
    o.name = name
    o.scale = radii
    bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)
    if material:
        o.data.materials.append(material)
    bpy.ops.object.shade_smooth()
    return o

def rod(name, p0, direction, length, r0, r1, material=None, verts=12):
    """Tapered rod from p0 along direction (rounded tip)."""
    d = Vector(direction).normalized()
    mid = Vector(p0) + d * (length / 2)
    bpy.ops.mesh.primitive_cone_add(vertices=verts, radius1=r0, radius2=r1, depth=length,
                                    location=mid)
    o = bpy.context.active_object
    o.name = name
    o.rotation_mode = 'QUATERNION'
    o.rotation_quaternion = Vector((0, 0, 1)).rotation_difference(d)
    bpy.ops.object.transform_apply(location=False, rotation=True, scale=False)
    if material:
        o.data.materials.append(material)
    bpy.ops.object.shade_smooth()
    return o

def smooth(vals, xs, sigma):
    out = []
    for xi in xs:
        w = [math.exp(-0.5 * ((xj - xi) / sigma) ** 2) for xj in xs]
        out.append(sum(wj * v for wj, v in zip(w, vals)) / sum(w))
    return out

def interp(ctrl, x):
    for (x0, *v0), (x1, *v1) in zip(ctrl, ctrl[1:]):
        if x0 <= x <= x1:
            t = (x - x0) / (x1 - x0)
            return [a + (b - a) * t for a, b in zip(v0, v1)]
    return ctrl[0][1:] if x < ctrl[0][0] else ctrl[-1][1:]

# ---------------- body (thorax + abdomen) ----------------
def build_body(P, m_body):
    L = P["seg_len"]
    x_abd0 = 1.00                       # thorax/abdomen junction
    a_abd = [0.45, 0.46, 0.46, 0.45, 0.43, 0.40, 0.36]   # half-widths I–VII
    b_abd = [0.31, 0.31, 0.30, 0.29, 0.28, 0.27, 0.25]   # half-heights I–VII
    ctrl = [(-0.12, 0.14, 0.12), (0.00, 0.24, 0.18), (0.15, 0.46, 0.31),
            (0.40, 0.60, 0.37), (0.68, 0.61, 0.37), (0.88, 0.54, 0.35), (x_abd0, 0.47, 0.32)]
    for i in range(7):
        ctrl.append((x_abd0 + L * (i + 0.5), a_abd[i], b_abd[i]))
    x8 = x_abd0 + 7 * L                  # start of segment VIII
    ctrl += [(x8 + 0.17, 0.31, 0.24), (x8 + 0.34, 0.22, 0.19),
             (x8 + 0.55, 0.16, 0.14), (x8 + 0.74, 0.10, 0.09), (x8 + 0.80, 0.04, 0.04)]
    x_end = x8 + 0.80

    n = 260
    xs = [ctrl[0][0] + (x_end - ctrl[0][0]) * i / (n - 1) for i in range(n)]
    A = [interp(ctrl, x)[0] for x in xs]
    B = [interp(ctrl, x)[1] for x in xs]
    A, B = smooth(A, xs, 0.05), smooth(B, xs, 0.05)
    bounds = [x_abd0] + [x_abd0 + L * i for i in range(1, 8)] + [x8 + 0.34]
    for k, x in enumerate(xs):
        f = 1.0
        for xb in bounds:
            f -= P["constriction"] * math.exp(-((x - xb) / 0.055) ** 2)
        A[k] *= f
        B[k] *= (1 - (1 - f) * 0.6)
    # taper the very last stations to close the tip smoothly
    A[-1], B[-1] = 0.02, 0.02

    N, e = P["ring_verts"], P["flatten_exp"]
    bm = bmesh.new()
    rings = []
    for x, a, b in zip(xs, A, B):
        ring = []
        for j in range(N):
            t = 2 * math.pi * j / N
            c, s = math.cos(t), math.sin(t)
            y = a * math.copysign(abs(c) ** (2 / e), c)
            z = b * math.copysign(abs(s) ** (2 / e), s)
            if z < 0:
                z *= 0.92          # ventral side slightly flatter
            ring.append(bm.verts.new((x, y, z)))
        rings.append(ring)
    for r0, r1 in zip(rings, rings[1:]):
        for j in range(N):
            bm.faces.new((r0[j], r0[(j + 1) % N], r1[(j + 1) % N], r1[j]))
    for ring, rev in ((rings[0], True), (rings[-1], False)):
        cx = sum((v.co for v in ring), Vector()) / N
        c = bm.verts.new(cx)
        for j in range(N):
            f = (ring[j], ring[(j + 1) % N], c)
            bm.faces.new(f[::-1] if rev else f)
    bmesh.ops.recalc_face_normals(bm, faces=bm.faces)
    me = bpy.data.meshes.new("body")
    bm.to_mesh(me)
    bm.free()
    o = bpy.data.objects.new("Body", me)
    bpy.context.collection.objects.link(o)
    o.data.materials.append(m_body)
    for p in o.data.polygons:
        p.use_smooth = True
    geo = dict(x_abd0=x_abd0, L=L, x8=x8, x_end=x_end, xs=xs, A=A, B=B)
    return o, geo

def half_size(geo, x):
    i = min(range(len(geo["xs"])), key=lambda k: abs(geo["xs"][k] - x))
    return geo["A"][i], geo["B"][i]

# ---------------- build ----------------
def build(P):
    for o in list(bpy.data.objects):
        bpy.data.objects.remove(o, do_unlink=True)
    m_body = mat("cuticle", (0.36, 0.30, 0.18), rough=0.35, trans=0.25)
    m_head = mat("head_capsule", (0.30, 0.21, 0.12), rough=0.3)
    m_dark = mat("eye_spiracle", (0.04, 0.03, 0.02), rough=0.25)
    m_app = mat("appendage", (0.42, 0.34, 0.22), rough=0.4, trans=0.3)
    t = P["appendage_thick"]
    parts = []

    body, g = build_body(P, m_body)
    parts.append(body)

    # head capsule: ovoid, slightly longer than wide, flattened
    hx = -0.33
    parts.append(ellipsoid("Head", (hx, 0, 0.02), (0.40, 0.31, 0.22), material=m_head, seg=32))
    # eyes (posterolateral)
    for s in (1, -1):
        parts.append(ellipsoid(f"Eye_{'L' if s > 0 else 'R'}", (hx + 0.17, s * 0.23, 0.07),
                               (0.07, 0.05, 0.06), material=m_dark, seg=16))
    # antennae: short, anterolateral, directed forward-outward
    for s in (1, -1):
        d = Vector((-math.cos(math.radians(30)), s * math.sin(math.radians(30)), 0.05))
        parts.append(rod(f"Antenna_{'L' if s > 0 else 'R'}", (hx - 0.27, s * 0.20, 0.06), d,
                         0.30, 0.035 * t, 0.018 * t, material=m_app))
    # lateral palatal (mouth) brushes: flattened fans, anteroventral
    for s in (1, -1):
        parts.append(ellipsoid(f"MouthBrush_{'L' if s > 0 else 'R'}",
                               (hx - 0.42, s * 0.10, -0.06),
                               (0.16, 0.08, 0.018 * t),
                               rot=(0, math.radians(-15), s * math.radians(25)),
                               material=m_app, seg=16))

    # spiracular plate on dorsum of VIII + two spiracles
    x8 = g["x8"]
    a8, b8 = half_size(g, x8 + 0.17)
    parts.append(ellipsoid("SpiracularPlate", (x8 + 0.17, 0, b8 * 0.88), (0.12, 0.15, 0.05),
                           material=m_head, seg=24))
    for s in (1, -1):
        parts.append(ellipsoid(f"Spiracle_{'L' if s > 0 else 'R'}", (x8 + 0.15, s * 0.05, b8 * 0.88 + 0.04),
                               (0.03, 0.025, 0.015), material=m_dark, seg=12))
    # saddle on segment X (dorsal sclerite, slightly raised)
    xX = x8 + 0.55
    aX, bX = half_size(g, xX)
    parts.append(ellipsoid("Saddle", (xX, 0, bX * 0.55), (0.13, aX * 0.95, bX * 0.55),
                           material=m_head, seg=20))
    # anal papillae: 4, posterior, splayed
    xe = g["x_end"] - 0.06
    for i, (sy, sz) in enumerate(((1, 1), (-1, 1), (1, -1), (-1, -1))):
        d = Vector((1, sy * 0.35, sz * 0.25)).normalized()
        c = Vector((xe, sy * 0.03, sz * 0.02)) + d * 0.15
        parts.append(ellipsoid(f"AnalPapilla_{i}", tuple(c), (0.16, 0.035 * t, 0.035 * t),
                               rot=Euler((0, -math.asin(d.z), math.atan2(d.y, d.x)), 'XYZ'),
                               material=m_app, seg=16))

    # simplified palmate hairs (seta 1) on III–VII
    if P["palmate"]:
        for i in range(2, 7):
            x = g["x_abd0"] + g["L"] * (i + 0.30)
            a, b = half_size(g, x)
            for s in (1, -1):
                parts.append(ellipsoid(f"Palmate_{i+1}_{'L' if s > 0 else 'R'}",
                                       (x, s * a * 0.55, b * 0.97),
                                       (0.05, 0.05, 0.012 * t), material=m_head, seg=12))
    # optional lateral setae
    if P["lateral_setae"]:
        for x, Ls in [(0.3, 0.55), (0.6, 0.6), (0.85, 0.5)] + \
                     [(g["x_abd0"] + g["L"] * (i + 0.4), 0.45 if i < 3 else 0.35) for i in range(7)]:
            a, b = half_size(g, x)
            for s in (1, -1):
                parts.append(rod("Seta", (x, s * a * 0.95, 0), (0.15, s, 0.05), Ls,
                                 0.008 * t, 0.003 * t, material=m_app, verts=6))
    return parts

def join(parts, name):
    bpy.ops.object.select_all(action='DESELECT')
    for o in parts:
        o.select_set(True)
    bpy.context.view_layer.objects.active = parts[0]
    bpy.ops.object.join()
    o = bpy.context.active_object
    o.name = name
    return o

def add_bend_modifiers(o):
    # origin at mid-body so bends are centered
    bpy.context.scene.cursor.location = (2.2, 0, 0)
    bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
    for nm, ax in (("Bend_lateral", 'Z'), ("Bend_dorsoventral", 'Y')):
        m = o.modifiers.new(nm, 'SIMPLE_DEFORM')
        m.deform_method = 'BEND'
        m.deform_axis = ax
        m.angle = 0.0

def watertight_copy(o, voxel, name):
    c = o.copy()
    c.data = o.data.copy()
    c.name = name
    bpy.context.collection.objects.link(c)
    bpy.ops.object.select_all(action='DESELECT')
    c.select_set(True)
    bpy.context.view_layer.objects.active = c
    for m in list(c.modifiers):
        bpy.ops.object.modifier_apply(modifier=m.name)
    r = c.modifiers.new("Remesh", 'REMESH')
    r.mode = 'VOXEL'
    r.voxel_size = voxel
    bpy.ops.object.modifier_apply(modifier="Remesh")
    d = c.modifiers.new("Dec", 'DECIMATE')
    d.ratio = 0.35
    bpy.ops.object.modifier_apply(modifier="Dec")
    bpy.ops.object.shade_smooth()
    return c

def export_stl(o, path, scale):
    bpy.ops.object.select_all(action='DESELECT')
    o.select_set(True)
    bpy.context.view_layer.objects.active = o
    bpy.ops.wm.stl_export(filepath=path, export_selected_objects=True, global_scale=scale,
                          apply_modifiers=True, ascii_format=False)

def export_obj(o, path):
    bpy.ops.object.select_all(action='DESELECT')
    o.select_set(True)
    bpy.context.view_layer.objects.active = o
    bpy.ops.wm.obj_export(filepath=path, export_selected_objects=True, export_materials=True,
                          apply_modifiers=True)

# ---------------- rendering ----------------
def setup_render(res=(1100, 600)):
    sc = bpy.context.scene
    sc.render.engine = 'CYCLES'
    sc.cycles.device = 'CPU'
    sc.cycles.samples = 48
    sc.render.resolution_x, sc.render.resolution_y = res
    sc.render.film_transparent = False
    w = sc.world or bpy.data.worlds.new("World")
    sc.world = w
    w.use_nodes = True
    w.node_tree.nodes["Background"].inputs["Color"].default_value = (0.92, 0.94, 0.95, 1)
    w.node_tree.nodes["Background"].inputs["Strength"].default_value = 0.5
    for loc, e in (((2, -4, 6), 25), ((-3, 4, 3), 8)):
        ld = bpy.data.lights.new("L", 'AREA')
        ld.energy = e
        ld.size = 4
        lo = bpy.data.objects.new("L", ld)
        lo.location = loc
        lo.rotation_euler = (Vector((2.2, 0, 0)) - Vector(loc)).to_track_quat('-Z', 'Y').to_euler()
        sc.collection.objects.link(lo)

def render_view(path, loc, target=(2.3, 0, 0), ortho=6.6):
    sc = bpy.context.scene
    cd = bpy.data.cameras.new("cam")
    cd.type = 'ORTHO'
    cd.ortho_scale = ortho
    co = bpy.data.objects.new("cam", cd)
    sc.collection.objects.link(co)
    co.location = loc
    up = 'Y'
    co.rotation_euler = (Vector(target) - Vector(loc)).to_track_quat('-Z', up).to_euler()
    sc.camera = co
    sc.render.filepath = path
    bpy.ops.render.render(write_still=True)
    bpy.data.objects.remove(co)

# ---------------- main ----------------
def main():
    argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(os.getcwd(), "out"))
    ap.add_argument("--print", action="store_true")
    ap.add_argument("--no-render", action="store_true")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    P = dict(PARAMS)
    if a.print:
        P["appendage_thick"] = 2.5
        P["palmate"] = False
    parts = build(P)
    larva = join(parts, "Anopheles_larva_L4")
    add_bend_modifiers(larva)
    if a.print:
        wt = watertight_copy(larva, P["voxel"], "print_mesh")
        export_stl(wt, os.path.join(a.out, f"anopheles_larva_print_x{int(P['print_scale'])}.stl"),
                   P["print_scale"])
        return
    wt = watertight_copy(larva, P["voxel"], "watertight_realsize")
    export_stl(wt, os.path.join(a.out, "anopheles_larva_realsize_mm.stl"), 1.0)
    wt.hide_render = True
    wt.hide_set(True)
    export_obj(larva, os.path.join(a.out, "anopheles_larva.obj"))
    if not a.no_render:
        setup_render()
        render_view(os.path.join(a.out, "render_dorsal.png"), (2.3, 0, 10))
        render_view(os.path.join(a.out, "render_lateral.png"), (2.3, -10, 0), ortho=6.6)
        render_view(os.path.join(a.out, "render_oblique.png"), (-1.5, -5.5, 4.5), ortho=6.8)
    bpy.ops.wm.save_as_mainfile(filepath=os.path.join(a.out, "anopheles_larva.blend"))

main()
