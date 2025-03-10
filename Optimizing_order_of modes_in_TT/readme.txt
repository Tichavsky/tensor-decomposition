Summary of the files for optimizing the order of modes
in a tensor train

dejkvadtt ..... given symmetric matrix Q and sampling range produces the sampled function as TT
numelall ...... computes number of parameters in a tensor train
prepoctiTT .... given tensor as a TT and new order of modes, computes new TT with reordered modes
ttsvd ......... svd-based computation of TT defined as an array of all its elements (Algorithm 1)
ttsvdtt.......  svd-based computation of TT given as a TT with non-ideal bond dimensions (Algorithm 2)
ttsvdopt .....  optimizing the order of modes of tensor defined as an array of all its elements (Algorithm 3)
ttcrossopt .... optimizing the order of modes of tensor given as TT, using TT-cross (Algorithm 4)
ttopttt .....  optimizing the order of modes of tensor given as TT, without using TT-cross




