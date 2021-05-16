'''
This class takes care of the routing when a website URL is selected based on the
flask framework

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
import re
import os
import json
import time
import string
import csv

from flask import render_template, jsonify, flash, redirect, url_for
from pyspark.sql import SparkSession
from pyspark.sql.functions import col, asc, desc
from scipy.spatial.distance import squareform
from sympy.core.tests.test_sympify import numpy

from app import app, data
from app.forms import SearchForm, SimilaritySearchForm, HeatmapForm
from app.ChemHammer_py import gen_numpy_arrays, min_flow_dist, min_flow_dist2
from app.MatrixSort import DistanceMatrixSorter

global global_csv
global_csv = ""

@app.route('/')
@app.route('/index')
def index():
    return render_template('index.html', title='Home')


@app.route('/atomic_search', methods=['GET', 'POST'])
def atomic_search():
    form = SearchForm()
    compounds = []
    scores = []
    time_taken = []
    total_count = 0

    if form.validate_on_submit():
        spark = SparkSession.builder \
            .appName('chemhammer_dev') \
            .getOrCreate()

        icsdFile = spark.read.parquet("processed_icsd")
        icsdFile.createOrReplaceTempView("icsdFile")
        icsdFile.printSchema()

        # rdd is used for brute-force search. Need to implement M-Tree
        # rdd还需要吗？csv文件怎么存储为tree？object为化学式？其他的属性为value?
        # 如果要搜索前50个相似的化合物，当第一个搜索完成以后，是重新开始搜索还是就近搜索？或者有某种顺序？

        rdd = icsdFile \
            .select('icsd_code', 'formula_structural', 'publication_date', 'journal_name', 'title', 'author',
                    'mod_petti_labels', 'ratios') \
            .filter( \
            col('formula_structural') \
                .contains(form.must_contain.data)) \
            .rdd \
            .repartition(2048)

        time_start = time.time()

        test_str = form.search_bar.data
        total_count = rdd.count()
        source_labels, source_demands = gen_numpy_arrays(test_str)

        scores = rdd.map(lambda x: (x["icsd_code"],
                                    x["formula_structural"],
                                    x["publication_date"],
                                    x["journal_name"],
                                    x["title"],
                                    x["author"],
                                    min_flow_dist(source_labels, source_demands, x['mod_petti_labels'], x['ratios']))) \
            .toDF(["Code", "Formula", "Date", "Journal", "Title", "Author", "Score"]) \
            .dropDuplicates() \
            .sort(col("Score") \
                  .asc()) \
            .head(int(form.results_to_display.data))

        time_taken = f"{time.time() - time_start:.2f}"

    return render_template('atomic_search.html', title='Search',
                           form=form,
                           compounds=compounds,
                           scores=scores,
                           total_count=total_count,
                           num_results=form.results_to_display.data,
                           time_taken=time_taken)


@app.route('/compound_similarity', methods=['GET', 'POST'])
def similarity_search():
    form = SimilaritySearchForm()
    distance = []
    time_taken = []

    if form.validate_on_submit():
        time_start = time.time()
        print(time_start)

        test_str1 = form.first_compound.data
        test_str2 = form.second_compound.data
        source_label1, source_demand1 = gen_numpy_arrays(test_str1)
        source_label2, source_demand2 = gen_numpy_arrays(test_str2)
        distance = min_flow_dist(source_label1, source_demand1, source_label2, source_demand2)

        time_taken = f"{time.time() - time_start:.2f}"
        print(time_taken)

    return render_template('compound_similarity.html',
                           form=form,
                           time_taken=time_taken,
                           distance=distance)


@app.route('/similarity_heatmap', methods=['GET', 'POST'])
def similarity_heatmap():
    form = HeatmapForm()
    distance_mat = []
    time_taken = []
    timestamp = []
    tsv_string = ""
    row_lb = []
    col_lb = []

    if form.validate_on_submit():
        time_start = time.time()
        compound1 = form.compound1.data
        compound2 = form.compound2.data
        compound3 = form.compound3.data
        compound4 = form.compound4.data
        compound5 = form.compound5.data
        compound6 = form.compound6.data
        compound7 = form.compound7.data
        compound8 = form.compound8.data
        compound9 = form.compound9.data
        compound10 = form.compound10.data
        formula_list = [compound1, compound2, compound3, compound4, compound5, compound6, compound7, compound8,
                        compound9, compound10]
        dis_array = []

        for i in range(len(formula_list[:-1])):
            for j in range(i + 1, len(formula_list)):
                dis_array.append(min_flow_dist2(formula_list[i], formula_list[j]))
                # print(dis_array)
        dis_arr = numpy.rint(dis_array)
        distance_mat = squareform(dis_arr)
        #print(distance_mat)

        sorted_distance_matrix = DistanceMatrixSorter(distance_mat).ordered_dist_mat
        print(sorted_distance_matrix)

        row_label = []
        column_label = []
        for i in range(10):
            for j in range(10):
                row_label.append(i)
                column_label.append(j)

        row_vals = []
        col_vals = []

        for i in range(10):
            row_lb.append(i+1)
            col_lb.append(i+1)
            for j in range(10):
                row_vals.append(i+1)
                col_vals.append(j+1)
                # row_lb.append(string.ascii_lowercase[i])
                # col_lb.append(string.ascii_uppercase[j])

        timestamp = time.time()

        distance_list = distance_mat.flatten()
        matrix = numpy.vstack((numpy.array(row_label).T, numpy.array(column_label).T, distance_list.T)).T



        list1 = zip(row_vals, col_vals, distance_list)

        tsv_string = "row_label,col_label,value\n"
        for row in list1:
            tsv_string += f"{row[0]},{row[1]},{row[2]}\n"

        data.csv_data = tsv_string
        print(len(tsv_string))

        global_csv = tsv_string
        print(global_csv)

        time_taken = f"{time.time() - time_start:.2f}"

    return render_template('similarity_heatmap.html', title='Similarity Heatmap',
                           form=form,
                           time_taken=time_taken,
                           distance_mat=distance_mat,
                           timestamp=timestamp,
                           csv=tsv_string,
                           row_lb=row_lb,
                           col_lb=col_lb
                           )


@app.route("/get-data", methods = ["GET", "POST"])
def returnData():
    f = data.csv_data
    return f
