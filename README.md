## MD analysis pipeline

Scripts required to process MD trajectories, to analyze peptide torsions, construct free energy landscapes and analyze them as described in **Molecular conformation as a determinant of tripeptide self-assembly**

## Environment variables

- Python 3.x
- MDAnalysis 2.7.0
- numpy 1.26.4
- matplotlib 3.9.1

## Example workflow

results/FDY1_C36m-sTIP3P/Run1
        -- example trajectory data for testing

python scripts/00_process_trajectory/ProcessTraj.py      
        -- correct trajectory PBC (example contains preprocessed trajectory)

python scripts/01_calculate_torsions/MeasureTorsions.py  
        -- measure torsional angles lambda1/lambda2

python scripts/02_make_FEL/RunFS.py                      
        -- construct FEL using 2D-KDE (run script invokes SolveFEL.py)

python scripts/03_analyze_FEL/AnalyzeFEL.py              
        -- evaluate FEL curvature 
