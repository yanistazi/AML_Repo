import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.model_selection import ShuffleSplit, GridSearchCV,train_test_split
from sklearn.model_selection import cross_validate
from sksurv.datasets import load_veterans_lung_cancer
from sksurv.column import encode_categorical
from sksurv.metrics import concordance_index_censored
from sksurv.svm import FastSurvivalSVM,FastKernelSurvivalSVM
from sksurv.kernels import clinical_kernel

def score_survival_model(model, X, y):
    prediction = model.predict(X)
    result = concordance_index_censored(y['Status'], y['Survival_in_days'], prediction)
    return result[0]
    
comp_ITD =[28]+list(range(189,205))

clin =list(range(158,165))
demo = [165,166]
comp_ITD_clin_demo = comp_ITD + clin + demo
all_gen_cyto = list(range(4,158))

df_final = pd.read_table("aml_prognosis_updated.tsv",sep="\t")
dict_features_type_final_comp = dict(zip(("all_gen_cyto"),
                                         (all_gen_cyto)))
estimator = FastSurvivalSVM(max_iter=1000, tol=1e-6, random_state=17)
param_grid = {'alpha': 10. ** np.array([-6,-5,-4,-3,-2,-1,0]),'optimizer':["avltree"]}
cv = ShuffleSplit(n_splits=5,random_state=17)
gcv = GridSearchCV(estimator, param_grid, scoring=score_survival_model,
                   n_jobs=50, iid=False, refit=True,
                   cv=cv)

df=pd.DataFrame(columns=dict_features_type_final_comp.keys())
#for key,item in dict_features_type_final_comp.items():
key = "comp_ITD_clin_demo"
item =comp_ITD_clin_demo
x = df_final.iloc[:,item]
y = np.array(list(zip(df_final.os_status, df_final.os)),dtype=[('Status', '?'), ('Survival_in_days', '<f8')])
ci=[]
for i in range(50):
    X_train, X_test, y_train, y_test = train_test_split(pd.DataFrame(x), y, test_size=0.2, random_state=i)
    gcv = gcv.fit(X_train,y_train)
    print(gcv.best_params_)
    ci.append(concordance_index_censored(y_test['Status'], y_test['Survival_in_days'], gcv.predict(X_test))[0])
    print(ci)
df[key] = ci
df.to_csv("SVM_comp_ITD_clin_demo.csv")    