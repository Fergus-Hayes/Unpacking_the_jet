# Unpacking the jet

A Bayesian analysis for constraining and comparing short gamma-ray jet structures given the gravitational wave and prompt emission observations of GW170817 and GW190425, as well at the rate of detected short gamma-ray bursts by the Swift detector as described in [Unpacking merger jets: a Bayesian analysis of GW170817, GW190425 and electromagnetic observations of short gamma-ray bursts](http://arxiv.org/abs/2305.06275). Top-hat (TH), Gaussian jet (GJ), power-law (PL) and double Gaussian (DG) jet structures are compared and the local rate of binary neutron star mergers is constrained.

## Python packages

An environment with both nessai and bibly are required. These can be installed via pip:

``` pip install nessai bilby ```

Or alternatively install the packages in `requirements.txt` via:

``` pip install -r requirements.txt ```

## Recreating results

To run the analysis for the Gaussian jet structure given GW170817 and a random seed of 73419243, and save to a local directory './outpath' use command:

``` python sample_model_nessai.py --path ./outpath --struct GJ --GW170817 --randomseed 73419243 ```

This run is represented in a tabulated form such that:

| Structure  | GW170817 | GW190425 | rates | fixed LF | random seed |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| GJ  | :heavy_check_mark: |  |  |  | 73419243  |

If instead we wish to run for the double Gaussian jet structure with GW170817, GW190425 and the rate of detected short gamma-ray bursts by Swift, and by using the fitted luminosity function, then run:

``` python sample_model_nessai.py --path ./outpath --struct DG --GW170817 --GW190425 --rates --Lfixed --randomseed 78951308 ```

This run is represented in the table as:

| Structure  | GW170817 | GW190425 | rates | fixed LF | random seed |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| DG  | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | 78951308  |

### Random seeds

The random seeds required to produce the results in the paper are given in the following table.

| Structure  | GW170817 | GW190425 | rates | fixed LF | random seed |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| TH | :heavy_check_mark: |  |  | :heavy_check_mark: | 81717499 |
| TH |  |  | :heavy_check_mark: | :heavy_check_mark: | 1898263 |
| TH |  | :heavy_check_mark: |  | :heavy_check_mark: | 81294950 |
| TH | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | 38433061 |
| TH | :heavy_check_mark: |  | :heavy_check_mark: | :heavy_check_mark: | 93999094 |
| TH | :heavy_check_mark: |  |  |  | 54520780 |
| TH |  |  | :heavy_check_mark: |  | 35185720 |
| TH |  | :heavy_check_mark: |  |  | 1898263 |
| TH | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |  | 8967996 |
| TH | :heavy_check_mark: |  | :heavy_check_mark: |  | 94081413 |
| GJ | :heavy_check_mark: |  |  | :heavy_check_mark: | 23096079 |
| GJ | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | 55748431 |
| GJ | :heavy_check_mark: |  | :heavy_check_mark: | :heavy_check_mark: | 46968036 |
| GJ |  |  | :heavy_check_mark: | :heavy_check_mark: | 12767752 |
| GJ |  | :heavy_check_mark: |  | :heavy_check_mark: | 73206677 |
| GJ | :heavy_check_mark: |  |  |  | 73419243 |
| GJ | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |  | 6231946 |
| GJ | :heavy_check_mark: |  | :heavy_check_mark: |  | 94081413 |
| GJ |  |  | :heavy_check_mark: |  | 35185720 |
| GJ |  | :heavy_check_mark: |  |  | 31277933 |
| PL |  | :heavy_check_mark: |  | :heavy_check_mark: | 669246 |
| PL | :heavy_check_mark: |  | :heavy_check_mark: | :heavy_check_mark: | 93999094 |
| PL | :heavy_check_mark: |  |  | :heavy_check_mark: | 69983429 |
| PL | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | 45040263 |
| PL |  |  | :heavy_check_mark: | :heavy_check_mark: | 70566340 |
| PL |  | :heavy_check_mark: |  |  | 1898263 |
| PL | :heavy_check_mark: |  | :heavy_check_mark: |  | 68158422 |
| PL | :heavy_check_mark: |  |  |  | 48197248 |
| PL | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |  | 91198413 |
| PL |  |  | :heavy_check_mark: |  | 35185720 |
| DG | :heavy_check_mark: |  | :heavy_check_mark: | :heavy_check_mark: | 23096079 |
| DG | :heavy_check_mark: |  |  | :heavy_check_mark: | 73206677 |
| DG | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | 78951308 |
| DG |  |  | :heavy_check_mark: | :heavy_check_mark: | 78408861 |
| DG |  | :heavy_check_mark: |  | :heavy_check_mark: | 79188414 |
| DG | :heavy_check_mark: |  | :heavy_check_mark: |  | 68158422 |
| DG | :heavy_check_mark: |  |  |  | 27815768 |
| DG | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |  | 29910191 |
| DG |  |  | :heavy_check_mark: |  | 35185720 |
| DG |  | :heavy_check_mark: |  |  | 12767752 |
