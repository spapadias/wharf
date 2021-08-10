import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import accuracy_score, f1_score

embeddings = pd.read_csv('model.w2v', header=None, delim_whitespace=True, skiprows=[0])
embeddings.set_index(0, inplace=True)

labels = pd.read_csv("data/labels/cora-labels", header=None, delim_whitespace=True)
labels.set_index(0, inplace=True)
labels.rename({1: 'label'}, axis=1, inplace=True)

dataset = embeddings.join(labels)
X = dataset.drop('label', axis=1)
y = dataset.label

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25)

clf = LogisticRegressionCV(Cs=10, cv=5, scoring="accuracy", verbose=False, max_iter=5000, solver='newton-cg')
clf.fit(X_train, y_train)
y_pred = clf.predict(X_test)

# print(f'Accuracy: {accuracy_score(y_test, y_pred)}')
# print(f'F1 Score (Macro): {f1_score(y_test, y_pred, average="macro")}')
# print(f'F1 Score (Micro): {f1_score(y_test, y_pred, average="micro")}')
# print(f'F1 Score (Weighted): {f1_score(y_test, y_pred, average="weighted")}')

with open("results.csv", "a") as myfile:
    myfile.write(f"{accuracy_score(y_test, y_pred)}, {f1_score(y_test, y_pred, average='macro')}, {f1_score(y_test, y_pred, average='micro')}, {f1_score(y_test, y_pred, average='weighted')} \n")