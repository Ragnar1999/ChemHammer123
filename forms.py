'''
Basic webforms to be rendered by the browser

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

from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField, BooleanField, SubmitField
from wtforms.validators import DataRequired


class SearchForm(FlaskForm):
    compound = StringField('Compound to Search', validators=[DataRequired()])
    results_to_display = StringField('Results to display', default='50', validators=[DataRequired()])
    must_contain = StringField('Must Contain Element (speeds up search when included)', default='')
    search = SubmitField('Search')


class SimilaritySearchForm(FlaskForm):
    first_compound = StringField('First Compound', validators=[DataRequired()])
    second_compound = StringField('Second Compound', validators=[DataRequired()])
    search = SubmitField('Calculate')


class HeatmapForm(FlaskForm):
    compound1 = StringField('Compound1', validators=[DataRequired()])
    compound2 = StringField('Compound2', validators=[DataRequired()])
    compound3 = StringField('Compound3', validators=[DataRequired()])
    compound4 = StringField('Compound4', validators=[DataRequired()])
    compound5 = StringField('Compound5', validators=[DataRequired()])
    compound6 = StringField('Compound6', validators=[DataRequired()])
    compound7 = StringField('Compound7', validators=[DataRequired()])
    compound8 = StringField('Compound8', validators=[DataRequired()])
    compound9 = StringField('Compound9', validators=[DataRequired()])
    compound10 = StringField('Compound10', validators=[DataRequired()])
    draw = SubmitField('Go')