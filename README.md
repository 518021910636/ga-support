# PGA / CGA Support Utilities

Small MATLAB examples of projective and conformal geometric algebra: basic
products, reversion, and elementary geometric representations.

This repository is a self-contained support-code sample for a technical
portfolio. It demonstrates algebra conventions and numerical checks. The
examples use synthetic coordinates and make no robotics performance claims.

## Included

- **PGA:** geometric product and reversion in the 3D Euclidean PGA even subalgebra.
- **CGA:** geometric / outer products, left contraction, grade selection,
  reversion, duality, and point / sphere / plane representations in Cl(4,1).
- One small example and deterministic tests, using base MATLAB only.

## Run

From this repository's root in MATLAB:

```matlab
addpath('matlab/pga', 'matlab/cga', 'examples', 'tests');
demo_support;
report = test_support;
```

The demo prints a squared Euclidean distance of `9` and a sphere-incidence
residual of `0`. Tests raise an error on failure and print a pass message.

## Conventions

`PGA_GeoProd(x,y)` returns an 8-by-1 column in the ordered basis
`[1, e23, e31, e12, e01, e02, e03, e0123]`, with `e0^2=0` and
`e1^2=e2^2=e3^2=1`. A six-component input represents only the six bivectors
and is embedded as `[0; x(:); 0]`. `PGA_Rev` expects eight components.
Reversion is not a general multivector inverse.

`cga5` uses 32 coefficients, indexed by blade bit mask plus one, in the
orthogonal basis `[e1,e2,e3,ep,em]` with metric `[1,1,1,1,-1]`.
It returns multivectors as rows. The null vectors are
`eo=(em-ep)/2` and `einf=ep+em`, so `eo.einf=-1`.
Points are `P(x)=eo+x+0.5*dot(x,x)*einf`.
A dual sphere is `P(c)-0.5*r^2*einf`; a dual plane is `n+d*einf`
for the equation `dot(n,x)=d`. Left contraction uses the convention
`<A_r B_s>_(s-r)` when `r<=s`, and zero otherwise.

Coordinates in the demo are dimensionless toy values. This small reference
implementation prioritizes readable conventions; no speed advantage is claimed.

## Scope and rights

The sample contains algebra support utilities and toy examples. Application
solvers, model data, and research results are outside this repository's scope.
See [RIGHTS.md](RIGHTS.md) for the retained file-specific license information.

---

本仓库用于作品集中的 PGA / CGA 基础代码展示：包含通用代数运算、
点/球/平面的表示、合成算例和测试。运行上述 MATLAB 命令即可验证。

作品集链接说明可写：**PGA / CGA 几何代数基础工具与数值验证示例（MATLAB）**。
