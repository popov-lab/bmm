# Construct a single finding for a pre-fit data report

Builds one finding in the shape
[`bmm_data_check()`](https://venpopov.com/bmm/dev/reference/bmm_data_check.md)
expects. Use it inside
[`data_check_findings()`](https://venpopov.com/bmm/dev/reference/data_check_findings.md)
methods rather than assembling the list by hand, so that an unrecognized
severity is caught where it is written instead of being silently printed
as a note.

## Usage

``` r
data_check_finding(severity, message)
```

## Arguments

- severity:

  Either `"warning"` for a likely data coding mistake, or `"note"` for
  something the user should merely be aware of. Warnings are listed
  before notes and printed in red.

- message:

  Character. The text shown in the report. It is re-wrapped when
  printed, so do not hard-wrap it.

## Value

A list with elements `severity` and `message`

## See also

[`data_check_findings()`](https://venpopov.com/bmm/dev/reference/data_check_findings.md),
[`bmm_data_check()`](https://venpopov.com/bmm/dev/reference/bmm_data_check.md)
