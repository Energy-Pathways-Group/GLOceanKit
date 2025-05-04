- Topic: Transformations

For constant stratification the vertical modes $$F$$ are cosines and thus computed with a discrete cosine transform (DCT). However, the vertical modes are normalized differently than the DCT.

These scaling factors are the inverse of the scaling factors from equations B12a and B12b in [Early, et al. (2021)](https://doi.org/10.1017/jfm.2020.995). The also include a sign normalization such that $F(0)>0$.

If `isHydrostatic==0` then the scaling factors do not include the $$f_0$$ and use the hydrostatic form of the eigendepths, `h`.
