# Security Closure Module

Evaluation metrics from the ISPD 2023 Contest.

## Commands

``` tcl
sec_evaluate
```

Evaluates the design for security metrics (placement and routing).

``` tcl
sec_read_design_metrics [-json filename]
```

- `-json`. Specify the JSON file to read.

Reads a JSON file containing the design metrics (power, setup_wns, hold_wns).

``` tcl
sec_read_baseline_metrics [-json filename]
```

- `-json`. Specify the JSON file to read.

Reads a JSON file containing the baseline metrics to compare the current design to.

``` tcl
sec_compute_overall
```

Computes the overall score based on ISPD23's metrics.

## FAQs

## License

BSD 3-Clause License. See [LICENSE](../../LICENSE) file.
