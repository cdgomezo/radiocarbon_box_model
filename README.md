# Time-series of GhG emissions estimation for the learning of inverse modeling

.. image:: https://zenodo.org/badge/351733532.svg
   :target: https://zenodo.org/badge/latestdoi/351733532

### Abstract

Inverse modeling is a commonly used method to infer greenhouse gases (GhGs) sources and sinks based on their observed concentrations in the atmosphere. It is a Bayesian framework and requires a priori fluxes of all the evaluated sources and sinks, atmospheric observations, an atmospheric transport model that relates the observations to the a priori fluxes and the uncertainties of both fluxes and observations. Various techniques exist to solve the inversion.

Atmospheric inverse modeling is and will become even more important in the future quantification of GhGs to monitor the compliance of the Nationally Determined Contributions (NDCs) under the Paris Agreement. Therefore, the scientific and educational communities are becoming more interested in using atmospheric inversions and this has risen a necessity of creating tools that facilitate understanding as well as training in these techniques.

Quantifying anthropogenic GhG emissions, such as $CO_2$ from fossil fuel burning or CH4 from human activities, from atmospheric concentration observations is difficult since the carbon from all sources, both natural and anthropogenic, is mixed in the atmosphere, making it necessary to use other signals or tracers to separate anthropogenic emissions from natural sources. For fossil fuel $CO_2$ emissions radiocarbon ($^{14}CO_2$) is an excellent emission tracer because, due to its radioactive decay (~ 5000 years), it cannot be found in fossil fuels, which have been deposited millions of years ago as organic material. We have developed a Jupyter Notebook based on Python for the quantification of multi-tracer GhGs fossil fuel emissions and its isotopes. The notebook solves for the emissions by applying atmospheric inversions within a practical two-box model. The inverse modeling notebook is based on the analytical maximum a posteriori (MAP) solution of the Bayes’ theorem and allows to assess the error in the state vector and its uncertainty.</p>

This basic but powerful notebook is meant to be an educational and training tool for university students and new researchers in the field as well as for researchers interested in the estimation of long-term (>centennial) time-series of GhG emissions since it is built as a modular algorithm to be easily modified, coupled or expanded to other approaches or models depending on the application. The notebook was initially developed for the inverse modeling of $CO_2$ and $^{14}CO_2$ simultaneously and it is being expanded for additional GHG such as $CH_4$ and $^{13}CH_4$.


