# Intro to Genomic Cross Prediction



Genomic prediction of the mean, variance and usefulness of crosses can be accomplished with `genomicMateSelectR` functions.

## Understanding mate selection

[**Click Here**](https://docs.google.com/presentation/d/1PW9OUOhk1ng1F2KXAO3FEq5sT1l_ioKYSiCEWY0rAkE/edit?usp=sharing) for a google slides presentation entitled "Genomic mate selection in outbred species: predicting cross usefulness with additive and total genetic covariance matrices" where I introduce the concepts.

The theory / formulae are summarized as part of [this `genomicMateSelectR` vignette](https://wolfemd.github.io/genomicMateSelectR/articles/non_additive_models.html#cross-predictions-with-genomicmateselectr-1)

### Literature

Here are some recommended articles to read regarding genomic mate selection (full citations will be at bottom): [@wolfe2021; @Bijma2020; @Bonk2016; @Lehermeier2017; @neyhart2019; @Neyhart2019; @werner2020]

In particular, read these:

1.  Wolfe et a. 2021.Genomic mating in outbred species: predicting cross usefulness with additive and total genetic covariance matrices. <https://doi.org/10.1093/genetics/iyab122>.
2.  Werner et al. 2020. Genomic selection strategies for clonally propagated crops. <https://doi.org/10.1101/2020.06.15.152017>.

You may very well need/want to read the literature that are referenced in these articles. If you do that, you'll have a solid foundation for understanding prediction of cross performance.

### Tutorial

For an a tutorial on how to execute these predictions using `genomicMateSelectR` functions, see this ["Getting started predicting crosses vignette"](https://wolfemd.github.io/genomicMateSelectR/articles/start_here.html).

## Non-additive effects

Up until now, we have used an additive-effects only model (`modelType="A"`), which gives us access to predictions of **GEBV**.

In addition, `genomicMateSelectR` enables two types of non-additive effects models to be implemented: an **additive plus dominance** model (`modelType="AD"`) and a directional dominance model that allows for an inbreeding depression (or heterotic) effect (`modelType="DirDom"`).

### Literature

-   For the basic quantitative genetics concepts of additive and dominance effects:

    -   [Intro to Quantitative Genetics](https://htmlpreview.github.io/?https://github.com/lfelipe-ferrao/lfelipe-ferrao.github.io/blob/master/class/survey/1.Introduction.html), part of Felipe Ferrão's "Survey of Breeding Tools and (Genomic Selection) Methods".
    -   See also the list of [Recommended Literature](recommended-literature) provided in the previous chapter [Intro to Genomic Prediction](intro-to-genomic-prediction)

-   See my summary of the additive and non-additive ["genomic prediction models implemented"](https://wolfemd.github.io/genomicMateSelectR/articles/non_additive_models.html#genetic-models-implemented-1) as part of the 2nd `genomicMateSelectR` vignette. See *also* the references to the literature, which are cited there.

-   Two good papers to start studying genomic prediction with non-additive effects are:

    -   Vitezica et al. 2013. "On the Additive and Dominant Variance and Covariance of Individuals Within the Genomic Selection Scope." *Genetics* 195 (4): 1223--30. <https://doi.org/10.1534/genetics.113.155176>

    -   Varona et al. 2013. "Non-Additive Effects in Genomic Selection." *Frontiers in Genetics* 9 (March). <https://doi.org/10.3389/fgene.2018.00078>

### Tutorial

The vignette in `genomicMateSelectR` entitled [**"Genomic prediction with non-additive effects"**](https://wolfemd.github.io/genomicMateSelectR/articles/non_additive_models.html)provides a complete tutorial on how to execute these models *and* predicting cross-performance with them.

## Parent-wise Cross-validation

**How can we estimate the accuracy of predicting previously untested crosses?**

In the mate selection article, @wolfe2021 I devised a cross-validation strategy that uses a pedigree -based approach, which I called "parent-wise cross-validation". The approach is described in detail in the manuscript. It is illustrated starting on [Slide 50 of *this* gSlides presentation](https://docs.google.com/presentation/d/1PW9OUOhk1ng1F2KXAO3FEq5sT1l_ioKYSiCEWY0rAkE/edit?usp=sharing).

`genomicMateSelectR` provides a function `runParentWiseCrossVal()` ([see here for the documentation / details](https://wolfemd.github.io/genomicMateSelectR/reference/runParentWiseCrossVal.html)) to implement this kind of cross-validation.

An example of it's implementation in-practice is part of the [IITA 2021 Genomic Selection documentation here](https://wolfemd.github.io/IITA_2021GS/05-CrossValidation.html#Parent-wise_cross-validation).

In the next section, I will attempt a smaller example using the data we have been working with in this manual.
