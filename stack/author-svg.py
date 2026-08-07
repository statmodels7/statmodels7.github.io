# Authors the hex stickers as hand-written SVG, one file per package, into
# logo/svg/. R keeps only the rasterization (make-logos.R). White chalk on
# the blackboard; the sanguine chalk survives only at the center of the
# statmodels7 honeycomb (Giovanni, 2026-08-07).
import math, os

INK, CHALK, EDGE = "#3D6B4C", "#F7F4D4", "#9C3E11"
SANG = "#DD7644"
FILLG = "#9BBBA2"
W, H, R = 520, 600, 250
CX, CY = W / 2, H / 2
NL = chr(10)

def hex_pts(cx, cy, r):
    pts = []
    for deg in range(90, 431, 60):
        a = math.radians(deg)
        pts.append((cx + r * math.cos(a), cy - r * math.sin(a)))
    return pts

def pline(pts, close=False):
    d = "M " + " L ".join(f"{x:.2f} {y:.2f}" for x, y in pts)
    return d + (" Z" if close else "")

def mapper(xlim, ylim, bx, by, bw, bh):
    def to(x, y):
        px = bx + (x - xlim[0]) / (xlim[1] - xlim[0]) * bw
        py = by + bh - (y - ylim[0]) / (ylim[1] - ylim[0]) * bh
        return px, py
    return to

def wordmark(name):
    size, avail = 46, 2 * (R * math.cos(math.radians(30))) - 52
    label = name + "7"
    while size > 18 and len(label) * 0.6 * size > avail:
        size -= 1
    return (f'<text x="{W//2}" y="408" font-family="Courier New,Courier,monospace" '
            f'font-size="{size}" font-weight="400" letter-spacing="0.5" '
            f'text-anchor="middle" fill="{CHALK}">{label}</text>' + NL)

def write(name, glyph):
    hexd = pline(hex_pts(CX, CY, R), close=True)
    svg = (f'<?xml version="1.0" encoding="UTF-8"?>' + NL +
           f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
           f'viewBox="0 0 {W} {H}">' + NL +
           f'<path d="{hexd}" fill="{INK}" stroke="{EDGE}" stroke-width="9"/>' + NL +
           glyph + wordmark(name) + "</svg>" + NL)
    os.makedirs("svg", exist_ok=True)
    with open(f"svg/{name}7.svg", "w", encoding="utf-8") as f:
        f.write(svg)
    print("wrote", name)

def stroke(d, col=CHALK, wpx=6.0, op=1.0, cap="round", dash=None):
    extra = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<path d="{d}" fill="none" stroke="{col}" stroke-width="{wpx}" '
            f'stroke-opacity="{op}" stroke-linecap="{cap}" '
            f'stroke-linejoin="round"{extra}/>' + NL)

def filled(d, col, op):
    return f'<path d="{d}" fill="{col}" fill-opacity="{op}" stroke="none"/>' + NL

def circle(x, y, r, fill, extra=""):
    return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{fill}"{extra}/>' + NL

def gamma_pdf(x, k, rate=1.0):
    if x <= 0: return 0.0
    return rate**k * x**(k-1) * math.exp(-rate*x) / math.gamma(k)

def norm_pdf(x, m, s):
    return math.exp(-0.5*((x-m)/s)**2) / (s*math.sqrt(2*math.pi))

# --- numericals7: a blackboard of the package's own symbols -----------------
def g_numericals():
    items = [
        # text, x, y, size, rotation, opacity
        ("∫", 200, 252, 124, -6, 1.0),     # integral
        ("Σ", 306, 206, 78, 5, 0.95),      # sum
        ("Γ", 362, 286, 56, -8, 0.8),      # gamma
        ("π", 158, 168, 52, -10, 0.7),     # pi
        ("∂", 286, 316, 54, -7, 0.7),      # partial
    ]
    g = ""
    for txt, x, y, size, rot, op in items:
        g += (f'<text x="{x}" y="{y}" font-family="Georgia,\'Times New Roman\',serif" '
              f'font-size="{size}" font-style="italic" text-anchor="middle" '
              f'fill="{CHALK}" fill-opacity="{op}" '
              f'transform="rotate({rot} {x} {y})">{txt}</text>' + NL)
    return g

# --- linkfunctions7: the logit alone, bold ----------------------------------
def g_linkfunctions():
    lg = lambda e: 1/(1+math.exp(-e))
    to = mapper((-5.2, 5.2), (0, 1), 124, 132, 272, 192)
    xs = [-5.2 + i*10.4/220 for i in range(221)]
    return stroke(pline([to(x, lg(x)) for x in xs]), CHALK, 9.0)

# --- distributions7: one density, filled ------------------------------------
def g_distributions():
    to = mapper((0, 11), (0, 0.23), 108, 122, 304, 204)
    xs = [i*11/260 for i in range(1, 261)]
    d1 = [(x, gamma_pdf(x, 4.0)) for x in xs]
    g = filled(pline([to(0.0, 0.0)] + [to(x, y) for x, y in d1] + [to(11, 0)],
                     close=True), FILLG, 0.75)
    g += stroke(pline([to(x, y) for x, y in d1]), CHALK, 7.0)
    g += stroke(pline([to(0, 0), to(11, 0)]), CHALK, 3.0, 0.8, cap="butt")
    return g

# --- optimizers7: descent into the bowl, all chalk --------------------------
def g_optimizers():
    cx, cy, ang = 260, 226, math.radians(-24)
    ca, sa = math.cos(ang), math.sin(ang)
    def E(a, b, n=140):
        pts = []
        for i in range(n+1):
            t = 2*math.pi*i/n
            x, y = a*math.cos(t), b*math.sin(t)
            pts.append((cx + x*ca - y*sa, cy + x*sa + y*ca))
        return pts
    g = ""
    for lev in range(5, 0, -1):
        a, b = 33*lev, 15*lev
        g += stroke(pline(E(a, b)), CHALK, 3.2, 0.42 + 0.10*(5-lev))
    lam = (1/33.0**2, 1/15.0**2)
    p = [148.0, 52.0]
    path = [tuple(p)]
    for _ in range(7):
        gv = (lam[0]*p[0], lam[1]*p[1])
        t = (gv[0]**2 + gv[1]**2) / (lam[0]*gv[0]**2 + lam[1]*gv[1]**2)
        p = [p[0]-t*gv[0], p[1]-t*gv[1]]
        path.append(tuple(p))
    pts = [(cx + x*ca - y*sa, cy + x*sa + y*ca) for x, y in path]
    g += stroke(pline(pts), CHALK, 4.5, 0.95, dash="10 7")
    for q in pts[:-1]:
        g += circle(q[0], q[1], 5.5, CHALK)
    g += circle(cx, cy, 8.0, CHALK)
    return g

# --- basis7: whole bells only, the central one filled ------------------------
def bspline_basis(knots, deg, i, x):
    if deg == 0:
        return 1.0 if knots[i] <= x < knots[i+1] else 0.0
    out = 0.0
    d1 = knots[i+deg] - knots[i]
    if d1 > 0: out += (x - knots[i]) / d1 * bspline_basis(knots, deg-1, i, x)
    d2 = knots[i+deg+1] - knots[i+1]
    if d2 > 0: out += (knots[i+deg+1] - x) / d2 * bspline_basis(knots, deg-1, i+1, x)
    return out

def g_basis():
    deg = 3
    knots = [k/8 for k in range(9)]
    to = mapper((0.0, 1.0), (0, 1.22), 95, 160, 330, 158)
    g = ""
    nb = len(knots) - deg - 1
    mid = nb // 2
    for i in range(nb):
        lo, hi = knots[i], knots[i+deg+1]
        xs = [lo + j*(hi-lo)/160 for j in range(161)]
        pts = [(x, bspline_basis(knots, deg, i, x)) for x in xs]
        if i == mid:
            g += filled(pline([to(lo, 0)] + [to(x, y) for x, y in pts] +
                              [to(hi, 0)], close=True), CHALK, 0.30)
            g += stroke(pline([to(x, y) for x, y in pts]), CHALK, 6.0)
        else:
            g += stroke(pline([to(x, y) for x, y in pts]), CHALK, 2.8, 0.65)
    g += stroke(pline([to(0.0, 0), to(1.0, 0)]), CHALK, 3.0, 0.8, cap="butt")
    return g

# --- parameters7: the rings of the chart, all chalk --------------------------
def g_parameters():
    A = (260, 128); B = (128, 330); C = (392, 330)
    def bary(p):
        return (p[0]*A[0] + p[1]*B[0] + p[2]*C[0],
                p[0]*A[1] + p[1]*B[1] + p[2]*C[1])
    def alr_inv(e1, e2):
        m = max(e1, e2, 0.0)
        u1, u2, u3 = math.exp(e1-m), math.exp(e2-m), math.exp(-m)
        s = u1+u2+u3
        return (u1/s, u2/s, u3/s)
    g = ""
    for r in [0.8, 1.7, 2.8, 4.2]:
        pts = [bary(alr_inv(r*math.cos(t), r*math.sin(t)))
               for t in [i*2*math.pi/180 for i in range(181)]]
        g += stroke(pline(pts, close=True), CHALK, 2.8, 0.75)
    cx0, cy0 = bary(alr_inv(0.0, 0.0))
    g += circle(cx0, cy0, 5.5, CHALK)
    g += stroke(pline([A, B, C], close=True), CHALK, 7.0)
    return g

# --- penalties7: three branches from one origin, weights tell them apart -----
def g_penalties():
    lam, a = 1.0, 3.7
    def scad(u):
        u = abs(u)
        if u <= lam: return lam*u
        if u <= a*lam: return (2*a*lam*u - u*u - lam*lam) / (2*(a-1))
        return lam*lam*(a+1)/2
    xmax = 5.4
    ymax = scad(xmax)*1.55
    to = mapper((-xmax, xmax), (0, ymax), 112, 128, 296, 192)
    xs = [-xmax + i*2*xmax/300 for i in range(301)]
    g = ""
    xr = [x for x in xs if 0.135*x*x <= ymax*0.97]
    g += stroke(pline([to(x, 0.135*x*x) for x in xr]), CHALK, 2.6, 0.5)
    g += stroke(pline([to(x, 0.62*abs(x)) for x in xs]), CHALK, 3.8, 0.7)
    g += stroke(pline([to(x, scad(x)) for x in xs]), CHALK, 7.5)
    g += stroke(pline([to(-xmax, 0), to(xmax, 0)]), CHALK, 3.0, 0.8, cap="butt")
    return g

# --- statmodels7: the honeycomb, the members named, the root at the center ---
def g_statmodels():
    # cells face each other vertex-to-vertex (they sit in the vertex
    # directions 30+60k), so the center distance must exceed 2*r
    r, cx, cy = 46, 260, 217
    dx = 2 * r + 10
    centers = [(cx, cy)]
    for k in range(6):
        a = math.radians(30 + 60*k)
        centers.append((cx + dx*math.cos(a), cy - dx*math.sin(a)))
    names = ["numericals", "linkfunctions", "distributions", "optimizers",
             "basis", "parameters", "penalties"]
    g = ""
    for i, (x, y) in enumerate(centers):
        cell = pline(hex_pts(x, y, r - 4), close=True)
        if i == 0:
            g += filled(cell, SANG, 0.9)
            g += stroke(cell, CHALK, 3.6)
        else:
            g += stroke(cell, CHALK, 3.6, 0.85)
        nm = names[i]
        size = 15
        while len(nm) * 0.62 * size > 2 * (r - 10) and size > 8:
            size -= 1
        col = INK if i == 0 else CHALK
        wgt = 700 if i == 0 else 400
        g += (f'<text x="{x:.1f}" y="{y + size*0.35:.1f}" '
              f'font-family="Courier New,Courier,monospace" font-size="{size}" '
              f'font-weight="{wgt}" text-anchor="middle" '
              f'fill="{col}">{nm}</text>' + NL)
    return g

# --- modelterms7: a scatter, the fit and its band, all chalk ------------------
def g_modelterms():
    to = mapper((0, 10), (-0.05, 1.15), 112, 126, 296, 200)
    def fit(x):
        return 0.30 + 0.34*math.sin(x*0.62 - 0.4) + 0.028*x
    def half(x):
        return 0.085 + 0.028*abs(x - 5.0)/5.0 + 0.04*(x/10)**2
    xs = [i*10/200 for i in range(201)]
    upper = [(x, fit(x) + half(x)) for x in xs]
    lower = [(x, fit(x) - half(x)) for x in reversed(xs)]
    g = filled(pline([to(x, y) for x, y in upper + lower], close=True),
               CHALK, 0.18)
    # the data: open chalk circles, deterministically jittered around the fit
    pts = []
    seedvals = [0.42, -0.83, 0.15, 0.91, -0.34, 0.67, -0.58, 0.23, -0.95,
                0.74, -0.12, 0.55, -0.71, 0.08, 0.88, -0.44, 0.31, -0.26]
    for k, u in enumerate(seedvals):
        x = 0.4 + k * 9.3 / (len(seedvals) - 1)
        y = fit(x) + u * half(x) * 1.35
        pts.append((x, y))
    for x, y in pts:
        px, py = to(x, y)
        g += circle(px, py, 5.5, "none",
                    extra=f' stroke="{CHALK}" stroke-width="2.6" stroke-opacity="0.9"')
    g += stroke(pline([to(x, fit(x)) for x in xs]), CHALK, 6.5)
    return g

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    write("numericals",    g_numericals())
    write("linkfunctions", g_linkfunctions())
    write("distributions", g_distributions())
    write("optimizers",    g_optimizers())
    write("basis",         g_basis())
    write("parameters",    g_parameters())
    write("penalties",     g_penalties())
    write("statmodels",    g_statmodels())
    write("modelterms",    g_modelterms())
