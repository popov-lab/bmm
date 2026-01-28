# DDM 7-Parameter Version Removal and 4-Parameter Stan Reversion - Summary

## Rationale

The 7-parameter version of the DDM has been removed and the Stan likelihood has been reverted to the 4-parameter version due to:

1. **Stan compatibility issues**: The 7par version requires the latest Stan version, but rstan is not keeping pace with Stan updates, creating compatibility problems.

2. **Computational feasibility**: The 7par version has infeasible run-times for typical use cases:
   - Single-subject models: 1-3 hours for 1500 trials
   - Hierarchical models: Days to weeks without HPC resources
   - Limited benefit from chain-level parallelization due to individual chain complexity

3. **Practical focus**: The 3par and 4par versions provide excellent performance and are sufficient for most research applications. The 4-parameter Stan likelihood (without variability parameters) provides a good balance between flexibility and computational efficiency.

## Changes Made

### Phase 1: Remove 7-Parameter Version

#### 1. Core Model Definition (`R/model_ddm.R`)
- **Removed**: Complete 7par entry from `.ddm_version_table` (original lines 73-109)
  - Removed parameters definition
  - Removed links configuration
  - Removed priors specification
  - Removed initialization ranges
- **Updated**: Default version in `.model_ddm()` function from "7par" to "4par"
- **Removed**: Legacy version alias `seven_par = "7par"` from version conversion

#### 2. Test Suite (`tests/testthat/test-model_ddm.R`)
- **Removed**: Basic construction test for 7par version
- **Removed**: Parameter definition test for 7par version
- **Removed**: Legacy version name test for "seven_par"
- **Removed**: Mock backend test for 7par version

#### 3. Documentation (`vignettes/articles/bmm_ddm.Rmd`)
- **Removed**: Computational warning about 7par infeasibility
- **Removed**: 7par runtime estimate from hierarchical model section
- **Removed**: Complete 7-parameter model section (~50 lines)
- **Removed**: 7par from model comparison examples

#### 4. News File (`NEWS.md`)
- **Updated**: Feature description now mentions only "3par" and "4par" versions

#### 5. Man Pages (`man/ddm.Rd`)
- **Auto-regenerated**: Using `devtools::document()`
- **Updated**: Default version now shows "4par" instead of "7par"

### Phase 2: Revert Stan Likelihood to 4-Parameter Version

#### 1. Stan Likelihood (`inst/stan_chunks/ddm_functions.stan`)
- **Changed**: Function signature from `ddm_lpdf(real rt, real mu, real drift, real bound, real ndt, real zr, real sdrift, real sndt, real szr, int dec)` to `ddm_lpdf(real rt, real mu, real drift, real bound, real ndt, real zr, int dec)`
- **Removed**: Variability parameters (sdrift, sndt, szr) from wiener_lpdf calls
- **Simplified**: Now calls `wiener_lpdf(rt | bound, ndt, zr, drift)` without variability parameters
- **Updated**: Comment changed from "7-parameter" to "4-parameter" ddm

#### 2. Family Definition (`R/model_ddm.R` - `configure_model.ddm`)
- **Removed**: link_sdrift, link_sndt, link_szr parameters from ddm_family function
- **Simplified**: dpars now only includes `c("mu","drift","bound","ndt","zr")` (removed sdrift, sndt, szr)
- **Simplified**: links now only includes 5 links (removed 3)
- **Simplified**: lb (lower bounds) now only has 5 elements (removed 3)
- **Simplified**: ub (upper bounds) now only has 5 elements (removed 3)
- **Removed**: link_sdrift, link_sndt, link_szr from formula$family construction

#### 3. Log-Likelihood Function (`R/model_ddm.R` - `log_lik_ddm`)
- **Removed**: szr, sdrift, sndt parameter extraction from prep object
- **Simplified**: Now only passes drift, bound, ndt, zr to dddm()

#### 4. Posterior Predict Function (`R/model_ddm.R` - `posterior_predict_ddm`)
- **Removed**: szr, sdrift, sndt parameter extraction from prep object
- **Simplified**: Now only passes drift, bound, ndt, zr to rddm()

#### 5. Model Definitions (`R/model_ddm.R` - `.ddm_version_table`)

**3par version**:
- **Removed**: sdrift, sndt, szr from parameters list (now only drift, bound, ndt, zr)
- **Removed**: sdrift, sndt, szr from links list
- **Simplified**: fixed_parameters now only contains `zr = 0.5` (removed sdrift = 0, sndt = 0, szr = 0)
- **Removed**: sdrift, sndt, szr from init_ranges

**4par version**:
- **Removed**: sdrift, sndt, szr from parameters list (now only drift, bound, ndt, zr)
- **Removed**: sdrift, sndt, szr from links list
- **Simplified**: fixed_parameters now empty list (no fixed parameters in 4par)
- **Removed**: sdrift, sndt, szr from init_ranges

#### 6. Test Updates (`tests/testthat/test-model_ddm.R`)
- **Updated**: 3par test to only check for zr fixed parameter (removed checks for sdrift, sndt, szr)
- **Updated**: 4par test to verify all 4 parameters exist and no fixed parameters (changed from checking variability parameters are fixed to 0)

## Verification

### Tests
- ✅ All DDM-specific tests passing (35 tests, reduced from 37 after removing 7par tests)
- ✅ No test failures related to changes
- ✅ No warnings in test suite

### R CMD CHECK
- ✅ 0 errors
- ✅ 0 warnings
- ✅ 3 benign notes (.claude directory, time verification, DDM_7PAR_REMOVAL_SUMMARY.md)

### Code Search
- ✅ Zero remaining references to "7par" in codebase
- ✅ All variability parameters (sdrift, sndt, szr) removed from Stan likelihood and family
- ✅ Model definitions updated to match Stan capabilities

## Current DDM Versions

After changes, the DDM model now supports:

### 3-Parameter Version (3par)
- **Parameters**: drift, bound, ndt
- **Fixed**: zr = 0.5 (symmetric starting point)
- **Stan likelihood**: Uses basic 4-parameter wiener_lpdf with fixed zr = 0.5
- **Use case**: Standard DDM with symmetric starting point
- **Performance**: Fast, suitable for large hierarchical models

### 4-Parameter Version (4par) - DEFAULT
- **Parameters**: drift, bound, ndt, zr
- **Fixed**: None (all parameters estimable)
- **Stan likelihood**: Uses 4-parameter wiener_lpdf
- **Use case**: DDM with flexible starting point bias
- **Performance**: Fast, suitable for hierarchical models
- **Link**: zr uses logit link (vs identity in 3par)

**Key change**: Both versions now use the same 4-parameter Stan likelihood function. The 3par version simply fixes zr = 0.5, while the 4par version estimates it. Neither version supports trial-to-trial variability parameters (sdrift, sndt, szr).

## Impact Assessment

### Breaking Changes
- Users calling `ddm(version = "7par")` will receive an error
- Users calling `ddm(version = "seven_par")` (legacy) will receive an error
- Any existing code using 7par must be updated to use 3par or 4par
- **Important**: Models fit before this change that used sdrift, sndt, or szr parameters will no longer work, as these parameters are no longer supported in the Stan likelihood

### Non-Breaking  
- Default version changed from 7par to 4par (internal function only)
- Public `ddm()` function already had 4par as default
- 3par and 4par functionality unchanged at the R interface level
- All existing 3par/4par formulas continue to work

### Stan Likelihood Simplification
- The Stan likelihood now only accepts 4 parameters: drift, bound, ndt, zr
- This matches the rtdists::ddiffusion implementation when variability parameters are set to 0
- R-side dddm() and rddm() functions still accept variability parameters for compatibility, but they default to 0
- When variability parameters are 0, the wiener_lpdf simplifies to the faster 4-parameter version

### Recommendations for Users
1. **For most applications**: Use 4par version (default)
2. **For symmetric designs**: Use 3par version for slightly faster computation
3. **For parameter recovery studies**: 3par or 4par provide good accuracy
4. **For trial-to-trial variability**: No longer supported - consider alternative models or implementations

## Files Modified

1. `inst/stan_chunks/ddm_functions.stan` - Simplified to 4-parameter likelihood
2. `R/model_ddm.R` - Model definitions, family, log_lik, and posterior_predict functions
3. `tests/testthat/test-model_ddm.R` - Updated test expectations
4. `vignettes/articles/bmm_ddm.Rmd` - Removed 7par documentation
5. `NEWS.md` - Updated feature description
6. `man/ddm.Rd` - Auto-regenerated documentation

## Technical Notes

- The Henrich et al. (2024) citation remains in the model documentation as it's relevant to the Stan implementation approach
- The R-side distribution functions (dddm, rddm) still accept variability parameters for API compatibility but default them to 0
- The Stan likelihood uses the simplified wiener_lpdf that doesn't include variability parameters
- Test coverage maintained with 35 tests (reduced from 37 after removing 7par-specific tests)
- All parameter definitions and fixed_parameters lists now match what's actually passed to Stan
