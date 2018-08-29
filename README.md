# Master code
This repository contains a set of Matlab programs for simulating superconducting phenomena in nanostructures. In particular, it is an implementation of the Usadel equation for diffusive nanostructures with superconducting, ferromagnetic, and spin-orbit coupled elements. The code was written as part of a Master project supervised by [Prof. Jacob Linder](https://folk.ntnu.no/jacobrun/) at NTNU, Norway. For more information about the assumptions and formalism underlying the code, please consult [the thesis itself](https://brage.bibsys.no/xmlui/handle/11250/2352094).

This code has been largely superseded by a [new program suite](https://github.com/jabirali/DoctorCode) written in modern object-oriented Fortran during my PhD. The new code is 2-3 orders of magnitude faster, can simulate a wider variety of physical systems, extends the software to handle non-equilibrium situations, and should also be easier to configure and use. Thus, the main reason for being interested in the code from my Master thesis, would be if you are a student that is planning to work on superconducting spintronics phenomena in Matlab, and could be interested in seeing how that can be done. Note that the implementation of strongly polarized spin-active interfaces was never rigorously tested and may therefore be prone to bugs.

The project itself is available under an MIT licence, which basically means that you can do anything you want with the code as long as you give credit where credit is due. The numerical solver in the subfolder `BVP/` was available for free under an unknown licence. Note that Matlab comes with another solver `bvp4c` which is largely compatible with the bundled `bvp6c`.
