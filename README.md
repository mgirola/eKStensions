# eKStensions

This project is based on the ROOT data analysis framework and was developed using * [ROOT CERN](https://root.cern/) v6.24/02 (Released 29th of June, 2021).
It should be compatible also with older versions as long as they where compiled using the C++17 standard.

The project consists in a self-explaining code meant to:
	1. Test one case in 1D where Kolmogorov-Smirnov test gives a result different from the chi^2;
	2. Test some of the provided Kolmogorov-Smirnov implementations that are currently found in ROOT, in particular the KolmogorovTest() attribute reimplemented in TH2* from TH1*;
	3. Implement and extend the Kolmogorov-Smirnov test from the 1D case to the 2D case in several ways, and perform some tests on them.
	4. Extend the Kolmogorov-Smirnov test to the case where the function has parameters which are estimated from the data themselves (this is not possible using a traditional approach, for more info search for: [Kolmogorov Smirnov](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test) and [Lilliefors Test](https://en.wikipedia.org/wiki/Lilliefors_test).

Contact me for any question.

