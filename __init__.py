'''
Initialised the components used 

Copyright (C) 2019  Cameron Hargreaves

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

--------------------------------------------------------------------------------
'''
import os

from flask import Flask, jsonify
from flask_bootstrap import Bootstrap
import json

app = Flask(__name__, static_url_path='/static')
app.config['SECRET_KEY'] = 'any secret string'
bootstrap = Bootstrap(app)


class DataStore():
    csv_data = None
    dist_mat_labels = None


data = DataStore()


from app import routes
