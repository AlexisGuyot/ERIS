# ERIS: Polarization, Antagonism and Community Boundaries

Documents in this repertory study the execution times of the matrix computation based version of ERIS on large graphs. A comparison is proposed with the previous version of ERIS (iterative) and Guerra et al.'s method (see https://github.com/rachel-bastos/boundaries-polarization).

-----------------

These documents include:
- **graphs/**: a repertory in which generated graphs are stored in GML format;
- **graph_rcis/**: a repertory in which the generated graphs for our article submitted to RCIS are stored in GML format;
- **times/**: a repertory in which the execution times measured are stored in CSV files;
- **times_rcis/**: a repertory in which the execution times measured for our article submitted to RCIS are stored in CSV files;
- **eris_iterative_complexity_analysis.R**: The R functions allowing to build the antagonism matrix with the iterative version of ERIS;
- **eris_iterative_complexity_analysis.R**: The R functions allowing to build the antagonism matrix with the matrix computation based version of ERIS;
- **guerra_complexity_analysis.py**: The Python functions allowing to build the antagonism matrix with Guerra et al.'s method;
- **Matrix ERIS Experiment - Execution on Large Graphs.ipynb**: A Jupyter Notebook explaining the methods used to generate the artificial graphs on which the methods were applied and to measured the execution times;
- **Matrix ERIS Experiment - Execution on Large Graphs - Python Part.ipynb**: A Jupyter Notebook measuring the execution times of the Python method;
- **Matrix ERIS Experiment - Execution on Large Graphs.html**: HTML version of the Jupyter Notebook;
- **Matrix ERIS Experiment - Execution on Large Graphs - Python Part.html**: HTML version of the Jupyter Notebook;
- **g_rt.gml**: The original anonymized graph built from real data harvested from Twitter and used to generate the artificial graphs.

-----------------

NB: For now, generated graphs for our article submitted to RCIS still need to be anonymized. Once this is done they will be made available in **graph_rcis/**.

Their caracteristics are:
| Vertices count | Edges count | Average out degree |
|----------------|-------------|--------------------|
| 684            | 681         | 0.99               |
| 1,669          | 1,686       | 1.01               |
| 2,654          | 2,710       | 1.02               |
| 3,912          | 4,127       | 1.05               |
| 20,604         | 29,390      | 1.42               |
| 35,516         | 74,217      | 2.08               |
| 53,780         | 106,975     | 1.99               |
| 72,054         | 159,127     | 2.21               |
| 106,896        | 229,541     | 2.15               |
| 180,981        | 433,222     | 2.39               |
| 373,997        | 1,164,409   | 3.11               |
| 673,887        | 2,535,423   | 3.76               |
| 904,458        | 4,427,838   | 4.90               |

-----------------

Plots of the execution times measured are available in **Matrix ERIS Experiment - Execution on Large Graphs.ipynb** and in its HTML counterpart.