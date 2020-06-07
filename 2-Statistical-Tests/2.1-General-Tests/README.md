# 2.1 Generic Tests

## Organization
---
- 2\.1\.1-Chi-Square-Test
    - Implementation of Knuth's offline Algorithm[1]
- 2\.1\.2-Kolmogorov-Smirnov-Test
    - Implementation of Sahni's O(n) offline Algorithm[2]
- 2\.1\.3-Piped
    - Piped implementation of 2.1.1 and 2.1.2
- 2\.1\.4-Quantile-Sketch
    - Online computation of quantiles using Greenwald's __Quantile Sketch__[3]
    - Implementation of Ashwin Lall's online algorithm for Kolmogorov-Smirnov Test[4]
    - My research report (unpublished) on __Streaming Algorithms for Chi-Square and Kolmogorov-Smirnov Tests using Quantile Sketch__

## References
---
1. Section 3.3.1, Donald E. Knuth, "The Art of Computer Programming Volume 2: Seminumerical Algorithms", Addison-Wesley Professional, 1997
2. Gonzalez, Teofilo F., Sartaj Sahni, and William R. Franta. "An efficient algorithm for the Kolmogorov-Smirnov and Lilliefors tests." ACM Trans. Math. Softw. 3.1 (1977): 60-64
3. Greenwald, Michael, and Sanjeev Khanna. "Space-efficient online computation of quantile summaries." ACM SIGMOD Record 30.2 (2001): 58-66
4. Lall, Ashwin. "Data streaming algorithms for the Kolmogorov-Smirnov test." 2015 IEEE International Conference on Big Data (Big Data). IEEE, 2015