# SimpleCentrifugalPump
Parametric design of centrifugal pump with single circular arc blade.

## Requirements

  * install the last stable version of `python3`


## Installation & Use

  * clone or download the code
  * open a terminal in `SimpleCentrifugalPump` and run:
    ```
    python3 scp.py [flow] [head]
    ```
    where `[flow]` rate in `[m^3/s]`  and `[head]` in `[m]` are the respective numeric value, for eg:
    ```
    python3 scp.py .011 25
    ```
  * it's also possible to specify other optional arguments, which they take a default value if not indicated, for eg:
    ```
    python3 scp.py .011 25 --hz=50 --t=.003
    ```
  * an help, in usage and in the meaning of the input arguments, can be obtaind with:
    ```
    python3 scp.py -h
    ```
