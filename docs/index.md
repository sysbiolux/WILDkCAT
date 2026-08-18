## Summary 

**[WILDkCAT](https://github.com/sysbiolux/WILDkCAT)** is a Python package that offers a comprehensive range of functions for extracting, retrieving and predicting enzyme kinetic data for a given metabolic model. The tool facilitates enzyme-constrained modelling by assigning specific kcat values to each combination of reaction, enzyme, and substrate(s).

WILDkCAT integrates the **[BRENDA](https://www.brenda-enzymes.org/)** and **[SABIO-RK](https://sabiork.h-its.org/)** APIs to retrieve experimental values and **[CataPro](https://github.com/zchwang/CataPro)** to predict missing values. It can be used with various enzyme-constrained metabolic modelling methods, such as **[COBREXA.jl](https://github.com/COBREXA/COBREXA.jl)** and **[GECKO](https://github.com/SysBioChalmers/GECKO)**, facilitating the seamless integration of kcat values into a given metabolic model.

---

<div style="text-align: justify">WILDkCAT also generates HTML reports at each step of the process, providing insights into the data retrieval and prediction outcomes: </div>

<p align="center"> <img src="report_example.gif" alt="WILDkCAT Report Demo" width="700"/> </p>

--- 

## Publication

WILDkCAT is described in our article published in _Bioinformatics_:

> Escoffier H, Linster CL, Sauter T. WILDkCAT: extract, retrieve, and predict enzyme turnover numbers of constraint-based metabolic models. Bioinformatics. (2026). [doi:10.1093/bioinformatics/btag510](doi:10.1093/bioinformatics/btag510)