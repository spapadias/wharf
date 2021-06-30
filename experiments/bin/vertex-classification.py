import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.metrics import accuracy_score

embeddings = pd.read_csv('model.w2v', header=None, delim_whitespace=True, skiprows=[0])
embeddings.set_index(0, inplace=True)

labels = pd.read_csv("data/labels/wiki-labels", header=None, delim_whitespace=True)
labels.set_index(0, inplace=True)
labels.rename({1: 'label'}, axis=1, inplace=True)

dataset = embeddings.join(labels)
X = dataset.drop('label', axis=1)
y = dataset.label

X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.2, test_size=None)

clf = LogisticRegression()
clf.fit(X_train, y_train)

y_pred = clf.predict(X_test)
acc = accuracy_score(y_test, y_pred)
print(acc)