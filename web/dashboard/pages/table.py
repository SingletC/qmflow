"""Instantiate a Dash app."""
import os
import time

import ase.db
import dash
from dash import dash_table, Output, State, Input
from dash import dcc
from dash import html
from rdkit import Chem
from rdkit.Chem import MolFromSmiles

from calc.call import SubmitTDDFTViaAndromeda
from web.dashboard.pages.data import create_dataframe
from web.dashboard.pages.layout import html_layout
import threading

from web.dashboard.pages.utils import split_filter_part
class CallBacks:
    def __init__(self, app: dash.Dash, db: ase.db.core.Database):
        self.db = db
        self.app = app
        self.df = create_dataframe(db)
        self.update_df_thead = threading.Thread(target=self.update_df)
        self.update_df_thead.start()

    def init_callbacks(self):

        # @self.app.callback(
        #     Output('database-table', 'data'),
        #
        #     Input('submit-val', 'n_clicks'),
        #
        # )
        # def update_table_data(n_clicks):
        #     """Create Dash datatable from Pandas DataFrame."""
        #     self.df = create_dataframe(self.db)
        #     return self.df .to_dict("records")

        @self.app.callback(
            Output('database-table', "data"),
            Input('database-table', "page_current"),
            Input('database-table', "page_size"),
            Input('database-table', "sort_by"),
            Input('database-table', "filter_query"), )
        def update_table(page_current, page_size, sort_by, filter):
            filtering_expressions = filter.split(' && ')
            dff = self.df.copy()
            for filter_part in filtering_expressions:
                if '{name}' in filter_part:
                    col_name = 'name'
                    if 'scontains' in filter_part:
                        operator = 'contains'
                        filter_value = filter_part.split(' ')[-1]
                    else:
                        col_name, operator, filter_value = split_filter_part(filter_part)
                else:
                    col_name, operator, filter_value = split_filter_part(filter_part)

                if col_name != 'name':
                    dff = dff.loc[getattr(dff[col_name], operator)(filter_value)]
                else:
                    if operator == 'contains':
                        filter_value = filter_value.replace("'", '')
                        match = Chem.MolFromSmiles(filter_value)
                        p = Chem.AdjustQueryParameters.NoAdjustments()
                        p.makeDummiesQueries = True
                        match = Chem.AdjustQueryProperties(match, p)
                        dff = dff[dff['name'].apply(lambda x: MolFromSmiles(x).HasSubstructMatch(match))]

                    elif operator == 'eq':
                        filter_value = filter_value.replace("'", '')
                        match = Chem.MolFromSmiles(filter_value)
                        canoical_smiles = Chem.MolToSmiles(match)
                        dff = dff[dff['name'] == canoical_smiles]

            if len(sort_by):
                dff = dff.sort_values(
                    [col['column_id'] for col in sort_by],
                    ascending=[
                        col['direction'] == 'asc'
                        for col in sort_by
                    ],
                    inplace=False
                )
            page = page_current
            size = page_size
            dff.ctime = dff.ctime.apply(lambda x: x.ctime())
            return dff.iloc[page * size: (page + 1) * size].to_dict('records')

        @self.app.callback(
            dash.dependencies.Output('database-table', 'page_current'),
            [dash.dependencies.Input('database-table', 'filter_query')])
        def reset_to_page_0(filter_query):

            return 0
    def update_df(self):
        while True:
            self.df = create_dataframe(self.db)
            time.sleep(120 if os.getenv('DEBUG') is not None else 10)


app = dash.register_page(__name__,path='/')
db = ase.db.connect(os.getenv('DATABASE_DIR'))
# Load DataFrame
df = create_dataframe(db)
df.ctime = df.ctime.apply(lambda x: x.ctime())
# Custom HTML layout
index_string = html_layout

# Create Layout
layout = html.Div(
    children=[
        # dcc.Graph(
        #     id="histogram-graph",
        #     figure={
        #         "data": [
        #             {
        #                 "x": df["complaint_type"],
        #                 "text": df["complaint_type"],
        #                 "customdata": df["key"],
        #                 "name": "311 Calls by region.",
        #                 "type": "histogram",
        #             }
        #         ],
        #         "layout": {
        #             "title": "NYC 311 Calls category.",
        #             "height": 500,
        #             "padding": 150,
        #         },
        #     },
        # ),

        dash_table.DataTable(
            id="database-table",
            columns=[
                {"name": i, "id": i, 'presentation': 'markdown' if i == 'Structure' else 'input', "deletable": True}
                for
                i in df.columns],
            data=df.to_dict("records"),
            sort_action='custom',
            sort_mode='multi',
            sort_by=[],
            page_size=20,
            page_current=0,
            page_action='custom',
            filter_action="custom",
            filter_query='',
            row_deletable=True,
            markdown_options={'link_target': '_blank', "html": True},
            style_cell_conditional=[{'if': {'column_id': 'Structure'},
                                     'width': '400px'},
                                    {'if': {'column_id': 'name'},
                                     'maxWidth': '50px'},
                                    ],

        ),
    ],
    id="dash-container",
)
callbacks = CallBacks(dash, db)
callbacks.init_callbacks()


