# Unpacking_the_jet

## Python packages

An environment with both nessai and bibly are required. These can be installed via pip:

``` pip install nessai bilby ```

## Recreating results

``` python sample_model_nessai.py --path ./outpath --struct GJ --GW170817 --randomseed 73419243 ```

| Structure  | GW170817 | GW190425 | rates | fixed LF | random seed |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| GJ  | :heavy_check_mark: |  |  |  | 73419243  |

``` python sample_model_nessai.py --path ./outpath --struct DG --GW170817 --GW190425 --rates --Lfixed --randomseed 78951308 ```

| Structure  | GW170817 | GW190425 | rates | fixed LF | random seed |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| DG  | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | 78951308  |

### Random seeds

| Structure  | GW170817 | GW190425 | rates | fixed LF | random seed |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| TH  | :heavy_check_mark: |  |  |  |  |
| TH  | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |  |   |
| TH  | :heavy_check_mark: |  | :heavy_check_mark: |  |   |
| TH  | | :heavy_check_mark: |  |  |   |
| TH  |  |  | :heavy_check_mark: |  |   |
| TH  | :heavy_check_mark: |  |  | :heavy_check_mark: |  |
| TH  | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |   |
| TH  | :heavy_check_mark: |  | :heavy_check_mark: | :heavy_check_mark: |   |
| TH  | | :heavy_check_mark: |  | :heavy_check_mark: |  |
| TH  |  |  | :heavy_check_mark: | :heavy_check_mark: |  |
| GJ  | :heavy_check_mark: |  |  |  | 73419243  |
| GJ  | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |  |   |
| GJ  | :heavy_check_mark: |  | :heavy_check_mark: |  |   |
| GJ  | | :heavy_check_mark: |  |  |   |
| GJ  |  |  | :heavy_check_mark: |  |   |
| GJ  | :heavy_check_mark: |  |  | :heavy_check_mark: |  |
| GJ  | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |   |
| GJ  | :heavy_check_mark: |  | :heavy_check_mark: | :heavy_check_mark: |   |
| GJ  | | :heavy_check_mark: |  | :heavy_check_mark: |  |
| GJ  |  |  | :heavy_check_mark: | :heavy_check_mark: |  |
| PL  | :heavy_check_mark: |  |  |  |  |
| PL  | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |  |   |
| PL  | :heavy_check_mark: |  | :heavy_check_mark: |  |   |
| PL  | | :heavy_check_mark: |  |  |   |
| PL  |  |  | :heavy_check_mark: |  |   |
| PL  | :heavy_check_mark: |  |  | :heavy_check_mark: |  |
| PL  | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |   |
| PL  | :heavy_check_mark: |  | :heavy_check_mark: | :heavy_check_mark: |   |
| PL  | | :heavy_check_mark: |  | :heavy_check_mark: |  |
| PL  |  |  | :heavy_check_mark: | :heavy_check_mark: |  |
| DG  | :heavy_check_mark: |  |  |  |  |
| DG  | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |  |   |
| DG  | :heavy_check_mark: |  | :heavy_check_mark: |  |   |
| DG  | | :heavy_check_mark: |  |  |   |
| DG  |  |  | :heavy_check_mark: |  |   |
| DG  | :heavy_check_mark: |  |  | :heavy_check_mark: |  |
| DG  | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | 78951308 |
| DG  | :heavy_check_mark: |  | :heavy_check_mark: | :heavy_check_mark: |   |
| DG  | | :heavy_check_mark: |  | :heavy_check_mark: |  |
| DG  |  |  | :heavy_check_mark: | :heavy_check_mark: |  |
