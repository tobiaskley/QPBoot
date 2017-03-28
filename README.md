# QPBoot
Model validation using Quantile Spectral Analysis and Parametric Bootstrap techniques

This pakage can be used for validating parametric time-series models. There is a demo available for checking if a GARCH(1,1)
model is suitable for DAX returns via 
```
demo("DAX")
```
The main method is ``qpBoot`` and for its ``model`` argument there are several predefined models (``getGARCH()``, or ``getARMA()`` for example).
