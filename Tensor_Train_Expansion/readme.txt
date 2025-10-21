List of software for tensor train expansion (for review in the IEEE Signal Processing Letters)

scalarprod ... scalar product of two tensors in the tensor train format
frobnorm ..... computes Frobenius norm of a tensor in the TT format
frobrozdil ... Frobenius norm of a difference of two tensors in TT format
plusTT ....... expression of a sum of two TTs as a single TT (increases the bond dimension)
hdotTT ....... expression of an elementwise (Hadamard) product of two TTs as a single TT (increases the bond dimension)

ttsvdtt ...   SVD-based TT decomposition in the TT format
ttsvdu  ...   SVD-based TT decomposition in the SOTT format
expandTT  ... computes expansion of a TT as a sum of TTs of bond dimension bmax

Variants of ttsvdtt (see the directory on Optimizing the order of tensor modes)

ttsvdtt0 ... version published in the IEEE SP Letters
ttsvdtt2 ... version submitted to icassp 2026
ttsvdtt3 ... version that can be faster when the bond dimensions are large

Variants of ttsvdu:

ttsvdu0 ...  version that applies plusTT followed by ttsvdtt3 (suitable for SOTT with a few elements)
ttsvdu ...   algorithm with linear complexity in the number of the elements in SOTT

