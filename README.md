# n1-simulator

Simulator for exploring N-of-1 studies.
Accompanies [FULL PAPER CITATION HERE], which contains a full description of the methods.


## Installation

The code is a single R file, which you can copy into your project and source:

```{r}
source('n1-simulator.R')
```

## Examples

See `examples/` in this repository.


## Function Reference

### n1_simulate

Runs a simulation of a single N-of-1 clinical trial.

#### Function Signature

```{r}
n1_simulate <- function(
    # REQUIRED ARGUMENTS
    n_treatments, n_blocks,
    baseline_initial,
    effect_size,
    tc_in, tc_out, tc_outcome,
    sd_baseline, sd_baseline, sd_obs,
    treatment_period, sampling_timestep, noise_timestep,
    
    # OPTIONAL ARGUMENTS
    treatment_mat_by_block = NULL,
    random_seed = NA,
    baseline_func = NULL,
    return_data_frame = TRUE
)
```

#### Return value



#### Arguments

##### `n_treatments`

<dl>
    <dt>Description</dt>
    <dd>
        The number of different treatments in the study.
    </dd>
    
    <dt>Type</dt>
    <dd>scalar integer</dd>
    
    <dt>Range</dt>
    <dd>`n_treatments >= 2`</dd>
</dl>

##### `n_blocks`

<dl>
    <dt>Description</dt>
    <dd>
        The number of blocks in the study, within which each treatment occurs once.
    </dd>

    <dt>Type</dt>
    <dd>scalar integer</dd>

    <dt>Range</dt>
    <dd>`n_blocks >= 1`</dd>
</dl>

##### `baseline_initial`

<dl>
    <dt>Description</dt>
    <dd>
        Initial value of baseline outcome random walk.
        Ignored if `!is.null(baseline_func)`.
    </dd>
    
    <dt>Type</dt>
    <dd>real scalar</dd>
    
    <dt>Range</dt>
    <dd>`is.finite(baseline_initial)`</dd>
</dl>

##### `effect_size`

<dl>
    <dt>Description</dt>
    <dd>
        Vector of effect sizes, one for each treatment.
    </dd>
    
    <dt>Type</dt>
    <dd>real vector</dd>
    
    <dt>Dimensions</dt>
    <dd>`length(effect_size) == n_treatments`</dd>
    
    <dt>Range</dt>
    <dd>`all(is.finite(effect_size))`</dd>
</dl>

##### `tc_in`

<dl>
    <dt>Description</dt>
    <dd>
        Vector of run-in time constants, one for each treatment.
    </dd>
    
    <dt>Type</dt>
    <dd>real vector</dd>
    
    <dt>Dimensions</dt>
    <dd>`length(tc_in) == n_treatments`</dd>
    
    <dt>Range</dt>
    <dd>`all(is.finite(tc_in))`</dd>
</dl>

##### `tc_out`

<dl>
    <dt>Description</dt>
    <dd>
        Vector of wash-out time constants, one for each treatment.
    </dd>
    
    <dt>Type</dt>
    <dd>real vector</dd>
    
    <dt>Dimensions</dt>
    <dd>`length(tc_out) == n_treatments`</dd>
    
    <dt>Range</dt>
    <dd>`all(is.finite(tc_out))`</dd>
</dl>

##### `sd_baseline`

<dl>
    <dt>Description</dt>
    <dd>
        Standard deviation of baseline drift process.
        After 1 unit of time, the amount of change in the baseline will be distributed with a standard deviation of this value.
    </dd>
    
    <dt>Type</dt>
    <dd>real scalar</dd>
    
    <dt>Range</dt>
    <dd>`sd_baseline >= 0`</dd>
</dl>
