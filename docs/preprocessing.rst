Python preprocessing
====================

The default backend is an independent NumPy/SciPy/Numba implementation. It does
not invoke the historical native kernels. Geometry and feature smoothing use
consistent FEM mass with weight 0.5. The curvature estimator uses a 1-mm signed-
distance grid, 1.5-mm derivative sigma, 6-mm band and gradient fallback threshold
0.2. Coordinates are not rounded. The tracking mesh is not replaced by the
smoothed curvature geometry.

The concave-positive minimum principal-curvature feature is
H - sqrt(max(H*H-K, 0)). The recommended robust-unit energy retains the robust
cortical Gaussian fit, with unweighted data terms, unit edge penalties and
regularization 1. Existing Python HJ-flux skeleton extraction and pruned-branch
cleanup are retained: endpoint flux percentile 5, gyral/sulcal pruning 25/15 mm
and unsupported reconstruction limit 5 mm2.

Between-gyri seeds default to the sulcal ROI intersection connected through
sulcal vertices to the two groups' shared interface. The simple intersection
remains an explicit option. Within a merged region, all sulcal intersection
vertices are eligible, with gyral-class vertices in the same union as targets.
Whole-hemisphere mode extends this policy to all valid cortical parcels and
applies the cortex mask when supplied. No named-sulcus atlas is used to select
an intraparietal basin automatically.

No-config workflows use PreprocessConfig.recommended(). The bare constructor
PreprocessConfig() preserves the historical explicit cotan baseline. CLI profile
overrides are explicit. All grid/mass/energy/HJ settings and input hashes are
recorded in pipeline metadata. Alternative independently implemented Python
methods remain available for explicit comparisons. Historical implementations
are kept on the legacy and dev-c++ branches; the Python package does not invoke
their executable kernels. The legacy-native backend and --tools-json option
have been removed. The legacy-real-even FOD convention identifier and existing
Python algorithm/seed-policy names are retained without changing their math.
