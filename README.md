# DIRECTEUR

<Introduction>  
We have proposed an in silico approach to comprehensively predict the functionalities of foods, encompassing even processed foods. This prediction is accomplished through the utilization of machine learning on biomedical big data. Our focus revolves around disease-related protein pathways, wherein we statistically evaluate how the constituent compounds collaboratively regulate these pathways.

#You can run the proposed method on jupyter notebook with the following procedures: 

#Aggregate component compound-protein interactions inferred by logistic regression analysis into food-target protein influence relationships.

- Run: food_protein.ipynb   
- Input: 
  - Compound-protein-score data format json  
  - Food-constituent compound-weight<sup>*1</sup>  data format json  
- Output: score/{result}<sup>*2</sup>

*1: Food-constituent compounds were collected from FooDB. Weights were assigned according to the number of foods containing the constituent compounds and were calculated using Equation 3 in “Revealing comprehensive food functionalities and mechanisms of action through machine learning”.  
*2: Food-target protein influence relationships score (Influence score) was calculated using Equation 4 in “Revealing comprehensive food functionalities and mechanisms of action through machine learning”.

#Convert the pre-established influential connections between foods and target proteins into functional relationships between foods and diseases. The target proteins designated by foods were aligned with the pathways linked to each disease and assessed for disease applicability using enrichment analysis. 

Run: Pathway_analysis.ipynb  
*pathway.py is a class; place it in the same folder as pathway_analysis.ipynb.
