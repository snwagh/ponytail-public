# [Ponytail](https://eprint.iacr.org/2020/451.pdf): Maliciously secure matrix multiplication supplementary code 

**Disclaimer:** This code is being released solely for the purpose of reducing duplication of academic work. The code contains most (only the online benchmarking was done using [MP-SPDZ](https://github.com/data61/MP-SPDZ)) of the components used in the [paper](https://eprint.iacr.org/2020/451.pdf). The authors will not be able to actively provide support or assistance in this but feel free to create pull requests if you have important improvements that might assist others.


### Table of Contents
- [Requirements](#requirements)
- [Source Code](#source-code)
    - [Building the code](#building)
    - [Running the code](#running)
- [References](#references)

---
#### Requirements
Most dependencies should be built from the `install.sh` script (other than basic ones which can be installed using `apt install *`)

---
#### Source Code
Most of the code resides in the src folder. The [Microsoft SEAL](https://github.com/microsoft/SEAL) backend was probably (sorry) version 3.3.0. Most of the relevant code is in `./src/`

---
#### Building the code
To build the executable `main`, run the following commands:
```
./install.sh
mkdir build; cd build
cmake ..
make -j
```
---
#### Running the code
```
./main 
```

### References
For most of the HE stuff, refer to the [Microsoft SEAL](https://github.com/microsoft/SEAL) library and the tutorials, they are extensive and very good. For technical aspects, the paper and the references there in should be a good resource. To cite this work: 
```
@inproceedings{chen2020maliciously,
  title={{Maliciously Secure Matrix Multiplication with Applications to Private Deep Learning}},
  author={Chen, Hao and Kim, Miran and Razenshteyn, Ilya and Rotaru, Drago\c{s} and Song, Yongsoo and Wagh, Sameer},
  booktitle={Advances in Cryptology -- ASIACRYPT},
  year={2020}
}
```

---
If there are any issues running the code, create a git issue but know that you're generally on your own. Good luck! 
