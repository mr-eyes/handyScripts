from plotly import __version__
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import plotly.graph_objs as go
import plotly.plotly as py
import pandas as pd
import sys
import os

if len(sys.argv) < 2:
    exit("python generate_boxes.py <newPairwise_tsv> <ortho/non 1/0>")

newPairwise_file = sys.argv[1]
matching = 1

if len(sys.argv) > 2:
    matching = int(sys.argv[2])

if matching:
    file_name = os.path.basename(newPairwise_file).split(".")[0] + "_ortho_boxes.html"
else:
    file_name = os.path.basename(newPairwise_file).split(".")[0] + "_non_ortho_boxes.html"


with open(newPairwise_file, "r") as pairwise:
    names = list(next(pairwise).strip().split()[4:])
    ln = len(names)

df = pd.read_csv(newPairwise_file, sep="\t")

boxPlots = []

print(df)
print("---------------------")
print("---------------------")

for i in range(0, ln, 1):
    #print("df " + str(i) +" done")
    _df = df.loc[df['L-' + str(i)] == matching]
    _sims = list(_df.iloc[0:, 0])
    print(_df)
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(_sims)
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    boxPlots.append(_sims)

traces = []
for i in range(len(boxPlots)):
    print("GO " + str(i) +" done")
    _name = "L-" + str(i+1)
    go_ob = go.Box(
        y=boxPlots[i],
        name=_name,
        marker=dict(
            color='rgb(8, 81, 156)',
        ),
        boxmean=True
    )
    traces.append(go_ob)


layout = go.Layout(
    yaxis=dict(
        type='log',
        autorange=True
    )
)

fig = dict(data=traces, layout=layout)
plot(fig,include_plotlyjs = True, filename= file_name, auto_open=True)

#py.iplot(traces, filename='WebGLmillion')
