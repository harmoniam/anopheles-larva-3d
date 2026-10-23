"""Orthographic three-view drawing (third-angle projection, A3, 40:1) of the larva model.
Requires: numpy, trimesh, scikit-image, matplotlib, and the Noto Sans CJK JP font.
Usage: python make_three_view.py [--obj path] [--out dir]
"""
import numpy as np, trimesh, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from skimage import measure
plt.rcParams["font.family"] = "Noto Sans CJK JP"

import argparse, os
ap = argparse.ArgumentParser(description="Make an A3 three-view drawing from the OBJ model.")
ap.add_argument("--obj", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "models", "anopheles_larva.obj"))
ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "docs"))
args = ap.parse_args()
sc = trimesh.load(args.obj, split_object=False, group_material=True)
V, F, M = [], [], []
mats = list(sc.geometry.keys())
off = 0
for mi, k in enumerate(mats):
    g = sc.geometry[k]
    v = g.vertices.copy()
    v = np.c_[v[:, 0], -v[:, 2], v[:, 1]]          # OBJ (Y-up) -> Blender (Z-up)
    V.append(v); F.append(g.faces + off); M.append(np.full(len(g.faces), mi)); off += len(v)
V, F, M = np.vstack(V), np.vstack(F), np.concatenate(M)

RES = 0.0025  # mm per pixel
def zbuffer(P, depth, ext):
    """P: (n,2) projected coords, depth: larger = closer to viewer."""
    (u0, v0), (u1, v1) = ext
    W, H = int((u1 - u0) / RES) + 1, int((v1 - v0) / RES) + 1
    Z = np.full((H, W), -np.inf); L = np.full((H, W), -1); S = np.zeros((H, W))
    tri = F
    pu = (P[:, 0] - u0) / RES; pv = (P[:, 1] - v0) / RES
    e1 = V[tri[:, 1]] - V[tri[:, 0]]; e2 = V[tri[:, 2]] - V[tri[:, 0]]
    nrm = np.cross(e1, e2); nrm /= np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-12
    return pu, pv, Z, L, S, nrm, W, H

def render(axis_u, axis_v, axis_d, light):
    P = np.c_[V @ axis_u, V @ axis_v]; D = V @ axis_d
    ext = ((P[:, 0].min() - .05, P[:, 1].min() - .05), (P[:, 0].max() + .05, P[:, 1].max() + .05))
    pu, pv, Z, L, S, nrm, W, H = zbuffer(P, D, ext)
    shade = np.clip(np.abs(nrm @ light), 0, 1)
    for t in range(len(F)):
        i = F[t]
        x, y, z = pu[i], pv[i], D[i]
        xa, xb = int(np.floor(x.min())), int(np.ceil(x.max())) + 1
        ya, yb = int(np.floor(y.min())), int(np.ceil(y.max())) + 1
        gx, gy = np.meshgrid(np.arange(xa, xb), np.arange(ya, yb))
        d = (y[1] - y[2]) * (x[0] - x[2]) + (x[2] - x[1]) * (y[0] - y[2])
        if abs(d) < 1e-12: continue
        w0 = ((y[1] - y[2]) * (gx - x[2]) + (x[2] - x[1]) * (gy - y[2])) / d
        w1 = ((y[2] - y[0]) * (gx - x[2]) + (x[0] - x[2]) * (gy - y[2])) / d
        w2 = 1 - w0 - w1
        ins = (w0 >= -1e-6) & (w1 >= -1e-6) & (w2 >= -1e-6)
        if not ins.any(): continue
        zz = w0 * z[0] + w1 * z[1] + w2 * z[2]
        gy_, gx_, zz = gy[ins], gx[ins], zz[ins]
        cur = Z[gy_, gx_]; upd = zz > cur
        Z[gy_[upd], gx_[upd]] = zz[upd]; L[gy_[upd], gx_[upd]] = M[t]; S[gy_[upd], gx_[upd]] = shade[t]
    return dict(Z=Z, L=L, S=S, ext=ext)

def contours(mask, ext):
    (u0, v0), _ = ext
    out = []
    for c in measure.find_contours(np.pad(mask.astype(float), 1), 0.5):
        c = c - 1
        out.append(np.c_[u0 + c[:, 1] * RES, v0 + c[:, 0] * RES])
    return out

ex, ey, ez = np.eye(3)
views = {
    "top":   render(ex, ey, ez, np.array([0.3, -0.3, 0.9])),       # 平面図 (dorsal)
    "front": render(ex, ez, -ey, np.array([0.3, -0.9, 0.4])),      # 正面図 (lateral, camera at -Y)
    "right": render(ey, ez, ex, np.array([0.9, 0.3, 0.4])),        # 右側面図 (from posterior)
}
from matplotlib.colors import ListedColormap
from skimage import measure
mats = list(mats)
RES = 0.0025; S = 40.0          # 40:1 on A3
CU = mats.index("cuticle")
# ---- measurements ----
top, front, right = views["top"], views["front"], views["right"]
def extent_mask(view, mask):
    (u0, v0), _ = view["ext"]
    ys, xs = np.nonzero(mask)
    return u0 + xs.min()*RES, u0 + xs.max()*RES, v0 + ys.min()*RES, v0 + ys.max()*RES
allx0, allx1, ally0, ally1 = extent_mask(top, top["L"] >= 0)
cux0, cux1, _, _ = extent_mask(top, top["L"] == CU)
def width_at(x, view=top, lab=None):
    (u0, v0), _ = view["ext"]
    col = int(round((x - u0) / RES))
    m = view["L"][:, col] >= 0 if lab is None else view["L"][:, col] == lab
    ys = np.nonzero(m)[0]
    return v0 + ys.min()*RES, v0 + ys.max()*RES
xs = np.arange(0.0, 1.0, 0.01); th = max(xs, key=lambda x: np.diff(width_at(x))[0]); thw = width_at(th)
x3 = 1.0 + 0.46*2.5; a3 = width_at(x3, lab=CU)
_, _, fz0, fz1 = extent_mask(front, front["L"] == CU)
_, _, rz0, rz1 = extent_mask(right, right["L"] >= 0)
ry0, ry1, _, _ = extent_mask(right, right["L"] >= 0)
head0, head1 = -0.73, 0.07
# print("len", allx1-allx0, "body", cux1-head0, "thorax w", thw, th, "abd III", a3, "h", fz0, fz1)

fig = plt.figure(figsize=(420/25.4, 297/25.4))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 420); ax.set_ylim(0, 297); ax.axis("off")
ax.add_patch(plt.Rectangle((10, 10), 400, 277, fill=False, lw=1.2))

OX = 35 - allx0*S            # paper X of model x=0
FY, TY, RX = 95, 190, 360    # front view z=0, top view y=0, right view y=0
cmap = ListedColormap(plt.cm.Greys(np.linspace(0.08, 0.45, 256)))

def place(view, fx, fy):
    (u0, v0), (u1, v1) = view["ext"]
    sh = np.where(view["L"] >= 0, 1 - view["S"], np.nan)
    ax.imshow(sh, origin="lower", cmap=cmap, vmin=0, vmax=1, interpolation="bilinear",
              extent=(fx(u0), fx(u1), fy(v0), fy(v1)), zorder=1)
    def draw(mask, lw, col="k"):
        for c in measure.find_contours(np.pad(mask.astype(float), 1), 0.5):
            c = c - 1
            ax.plot(fx(u0 + c[:, 1]*RES), fy(v0 + c[:, 0]*RES), color=col, lw=lw, zorder=3)
    for i in range(len(mats)):
        if i != CU: draw(view["L"] == i, 0.35)
    draw(view["L"] >= 0, 0.9)

fxT = lambda u: OX + u*S
place(top,   fxT, lambda v: TY + v*S)
place(front, fxT, lambda v: FY + v*S)
place(right, lambda u: RX + u*S, lambda v: FY + v*S)

# center lines
cl = dict(color="k", lw=0.35, ls=(0, (12, 3, 2, 3)), zorder=2)
ax.plot([OX + allx0*S - 6, OX + allx1*S + 6], [TY, TY], **cl)
ax.plot([RX, RX], [FY + rz0*S - 6, FY + rz1*S + 6], **cl)
# projection lines (thin)
pl = dict(color="0.6", lw=0.25, ls=":", zorder=0)
for x in (allx0, allx1):
    ax.plot([OX + x*S]*2, [FY + fz1*S + 3, TY + ally0*S - 3], **pl)
for z in (rz0, rz1):
    ax.plot([OX + allx1*S + 3, RX + ry0*S - 3], [FY + z*S]*2, **pl)

def dim(p1, p2, off, text, vertical=False, fs=8):
    (x1, y1), (x2, y2) = p1, p2
    kw = dict(color="k", lw=0.35, zorder=4)
    if not vertical:
        y = max(y1, y2) + off if off > 0 else min(y1, y2) + off
        ax.plot([x1, x1], [y1, y + np.sign(off)*1.5], **kw); ax.plot([x2, x2], [y2, y + np.sign(off)*1.5], **kw)
        ax.annotate("", (x1, y), (x2, y), arrowprops=dict(arrowstyle="<|-|>", lw=0.4, mutation_scale=6, shrinkA=0, shrinkB=0))
        ax.text((x1+x2)/2, y + 1.0, text, ha="center", va="bottom", fontsize=fs)
    else:
        x = max(x1, x2) + off if off > 0 else min(x1, x2) + off
        ax.plot([x1, x + np.sign(off)*1.5], [y1, y1], **kw); ax.plot([x2, x + np.sign(off)*1.5], [y2, y2], **kw)
        ax.annotate("", (x, y1), (x, y2), arrowprops=dict(arrowstyle="<|-|>", lw=0.4, mutation_scale=6, shrinkA=0, shrinkB=0))
        ax.text(x + (1.2 if off > 0 else -1.2), (y1+y2)/2, text, ha="left" if off > 0 else "right",
                va="center", fontsize=fs, rotation=90 if abs(y2-y1) < 30 and False else 0)

P = lambda x, y, Y0: (OX + x*S, Y0 + y*S)
# top view dims
dim(P(allx0, ally1, TY), P(allx1, ally1, TY), 22, f"全長 {allx1-allx0:.2f}（触角〜肛門乳頭）")
dim(P(head0, ally1, TY), P(cux1, ally1, TY), 12, f"体長 {cux1-head0:.2f}（頭部前縁〜X節末端）")
dim(P(head0, 0.31, TY), P(head1, 0.31, TY), 5, f"頭長 {head1-head0:.2f}", fs=7)
dim(P(0.0, thw[1], TY), P(1.0, thw[1], TY), 5, "胸部 1.00", fs=7)
L = 0.46
dim(P(1.0+L*3, a3[1], TY), P(1.0+L*4, a3[1], TY), 5, "IV節 0.46", fs=7)
dim(P(th, thw[0], TY), P(th, thw[1], TY), -0, "", vertical=True)
ax.texts[-1].remove() if False else None
# vertical width dims left of top view
def vdim(x, y0, y1, Y0, xpaper, text):
    kw = dict(color="k", lw=0.35, zorder=4)
    for y in (y0, y1):
        ax.plot([OX + x*S, xpaper - 1.5], [Y0 + y*S]*2, **kw)
    ax.annotate("", (xpaper, Y0 + y0*S), (xpaper, Y0 + y1*S),
                arrowprops=dict(arrowstyle="<|-|>", lw=0.4, mutation_scale=6, shrinkA=0, shrinkB=0))
    ax.text(xpaper + 1.2, Y0 + (y0+y1)/2*S, text, ha="left", va="center", fontsize=8, rotation=90)
ax.lines.pop() if False else None
vdim(th, thw[0], thw[1], TY, OX + allx0*S - 10, f"胸部最大幅 {thw[1]-thw[0]:.2f}")
ax.texts[-1].set_ha("right"); ax.texts[-1].set_x(OX + allx0*S - 11.5)
vdim(x3, a3[0], a3[1], TY, OX + x3*S, "")
ax.text(OX + x3*S + 1.5, TY + a3[0]*S - 4, f"III節幅 {a3[1]-a3[0]:.2f}", fontsize=7, ha="left", va="top")
# segment labels
segs = [("I",1.23),("II",1.69),("III",2.15),("IV",2.61),("V",3.07),("VI",3.53),("VII",3.99),("VIII",4.39),("X",4.78)]
for s, x in segs:
    ax.text(OX + x*S, TY + ally0*S - 4, s, ha="center", va="top", fontsize=7.5)
ax.text(OX + (-0.33)*S, TY + ally0*S - 4, "頭部", ha="center", va="top", fontsize=7.5)
ax.text(OX + 0.5*S, TY + ally0*S - 4, "胸部", ha="center", va="top", fontsize=7.5)
# front view dims
vdim(1.0+L*0.5, fz0, fz1, FY, OX + allx0*S - 10, f"最大高 {fz1-fz0:.2f}")
ax.texts[-1].set_ha("right"); ax.texts[-1].set_x(OX + allx0*S - 11.5)
# right view dims
dim((RX + ry0*S, FY + rz0*S), (RX + ry1*S, FY + rz0*S), -9, "")
ax.text(RX, FY + rz0*S - 13, f"最大幅 {ry1-ry0:.2f}", ha="center", va="top", fontsize=8)
kw = dict(color="k", lw=0.35, zorder=4)
xr = RX + ry1*S + 8
for z in (rz0, rz1): ax.plot([RX + 0.2*S, xr + 1.5], [FY + z*S]*2, **kw)
ax.annotate("", (xr, FY + rz0*S), (xr, FY + rz1*S), arrowprops=dict(arrowstyle="<|-|>", lw=0.4, mutation_scale=6, shrinkA=0, shrinkB=0))
ax.text(xr + 1.5, FY + (rz0+rz1)/2*S, f"最大高 {rz1-rz0:.2f}", rotation=90, va="center", ha="left", fontsize=8)
# part callouts (front/top)
def call(xy, txt, dxy):
    ax.annotate(txt, xy=xy, xytext=(xy[0]+dxy[0], xy[1]+dxy[1]), fontsize=7,
                arrowprops=dict(arrowstyle="-", lw=0.3), ha="center", va="center", zorder=5)
call(P(-0.62, 0.0, FY), "触角", (-4, 18))
call(P(-0.55, -0.12, FY), "口刷毛", (2, -20))
call(P(-0.16, 0.13, FY), "複眼", (6, 20))
call(P(4.39, 0.24, FY), "気門板", (0, 16))
call(P(4.78, 0.10, FY), "saddle", (8, 18))
call(P(5.1, -0.05, FY), "肛門乳頭", (6, -20))
call(P(2.3, 0.30, FY), "palmate hair（簡略）", (0, 20))

# view titles
ax.text(OX + 2.2*S, TY + ally1*S + 34, "平面図（背面）", ha="center", fontsize=11)
ax.text(OX + 2.2*S, FY + fz0*S - 32, "正面図（側面）", ha="center", fontsize=11)
ax.text(RX, FY + rz1*S + 22, "右側面図（後端から）", ha="center", fontsize=11)
vdim2 = None

# scale bar
sbx, sby = 40, 32
ax.plot([sbx, sbx + S], [sby, sby], color="k", lw=2.2, solid_capstyle="butt")
ax.text(sbx + S/2, sby + 2.5, "1 mm（実寸）", ha="center", fontsize=8)

# title block
tbx, tby, tbw, tbh = 260, 10, 150, 48
ax.add_patch(plt.Rectangle((tbx, tby), tbw, tbh, fill=False, lw=0.9))
rows = [("図名", "Anopheles 4齢幼虫 3Dモデル 三面図"), ("投影法", "第三角法"),
        ("尺度", "40 : 1（A3・100%印刷時）"), ("単位", "mm（実寸換算）"),
        ("注記", "剛毛は省略、palmate hair は円盤で簡略表示"),
        ("モデル", "anopheles_larva.blend")]
h = tbh / len(rows)
for i, (k, v) in enumerate(rows):
    y = tby + tbh - (i+1)*h
    if i: ax.plot([tbx, tbx+tbw], [y+h, y+h], color="k", lw=0.35)
    ax.text(tbx + 2, y + h/2, k, va="center", fontsize=7.5)
    ax.text(tbx + 20, y + h/2, v, va="center", fontsize=8)
ax.plot([tbx+18, tbx+18], [tby, tby+tbh], color="k", lw=0.35)

fig.savefig(os.path.join(args.out, "anopheles_larva_three_view.pdf"))
fig.savefig(os.path.join(args.out, "anopheles_larva_three_view.png"), dpi=200)
print("saved to", args.out)
