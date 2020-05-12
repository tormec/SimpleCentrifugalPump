# SimpleCentrifugalPump
Centrifugal pump parametric design with simple arc blades.

## Requirements

  * install the last stable version of `python3`


## Installation & Use

  * clone or download the code
  * open a terminal in `SimpleCentrifugalPump` and run:
    ```
    python3 scp.py [flowrate] [head]
    ```
    where `[flowrate]` in `[m^3/s]`  and `[head]` in `[m]` are the respective numeric value, for eg:
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
