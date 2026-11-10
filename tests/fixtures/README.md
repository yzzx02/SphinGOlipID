# Synthetic mzML fixture

`untargeted_two_local_peaks.mzML` is generated test data, not an experimental
sample: 161 centroided positive MS1 scans, two Gaussian signals of the same
measured m/z and three acquired MS2 scans. The first two MS2 have deliberately
simple Cer fragments. The third lacks required Cer evidence.

Default XCMS settings detect the first chromatographic peak only. This tests
that the second acquisition does not create a replacement feature when local
MS1 detection is absent. Separate unit fixtures test both independently
detected local peaks. No benchmark annotations seed this fixture or the library.
