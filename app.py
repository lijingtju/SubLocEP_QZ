from flask import Flask
from flask import render_template, url_for, flash, redirect, request, abort, send_file
import os
from pyecharts import Line
import pandas as pd
from forms import DataForm
app = Flask(__name__)
import os
import EP3_1
import EP3_2
import datetime
import os
nowTime = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')

SECRET_KEY = os.urandom(32)
app.config['SECRET_KEY'] = SECRET_KEY
REMOTE_HOST = "https://pyecharts.github.io/assets/js"

@app.route('/SubLocEP/WebServer',methods=['GET','POST'])
def hello1():
    form = DataForm()

    if request.method == 'POST':
        print(form.EP.data)
        f = request.files['dataset']
        import uuid
        uuid = uuid.uuid4().hex[-6:]
        form.dataset.data.filename = uuid + form.dataset.data.filename
        form.dataset.data.filename = str(nowTime) + "_" + form.dataset.data.filename
        pathname = f'/home/lijing/SubLocEP/Data/{form.dataset.data.filename}'
        f.save(pathname)
        print(pathname)
        if form.EP.data == "EP1":
            EP3_1.main(f'{form.dataset.data.filename}')

            # result_name=f'{pathname}_results.txt'
            return render_template('result_download.html', fileurl=form.dataset.data.filename + '_results.csv')
        EP3_2.main(f'{form.dataset.data.filename}')
        result_name = f'{pathname}_results.txt'
        return render_template('result_download.html', fileurl=form.dataset.data.filename + '_results.csv')
    return render_template('MRMD2.html', title='WebServer', form=form)




@app.route('/SubLocEP')
@app.route('/SubLocEP/Home')
def Home():
    return render_template('Home.html',title='Home')

@app.route('/SubLocEP/Download')
def Download():
    return render_template('Download.html',title='Download')

@app.route('/SubLocEP/Dataset')
def Dataset():
    return render_template('Dataset.html',title='Dataset')

@app.route('/SubLocEP/Help')
def Help():
    return render_template('Help.html',title='Help')

@app.route('/SubLocEP/Contact')
def Contact():
    return render_template('Contact.html', title='Contact')


@app.route("/SubLocEP/result_example")
def result_example():
    features_kv = [('uncharger.uncharger.gap3', '0.0095'),
                   ('negativecharger.alphaticr.gap3', '0.0093'),
                   ('aromatic.uncharger.gap4', '0.0092'),
                   ('postivecharger.alphaticr.gap3', '0.0090'),
                   ('alphaticr.alphaticr.gap4', '0.0089'),
                   ('alphaticr.negativecharger.gap3', '0.0089'),
                   ('postivecharger.uncharger.gap4', '0.0089'),
                   ('postivecharger.postivecharger.gap4', '0.0087'),
                   ('aromatic.alphaticr.gap4', '0.0086'), ('uncharger.aromatic.gap3', '0.0085'),
                   ('negativecharger.alphaticr.gap4', '0.0083'),
                   ('postivecharger.aromatic.gap3', '0.0082'),
                   ('postivecharger.negativecharger.gap5', '0.0082'),
                   ('alphaticr.negativecharger.gap5', '0.0082'),
                   ('aromatic.postivecharger.gap4', '0.0081'),
                   ('uncharger.aromatic.gap5', '0.0081'), ('aromatic.uncharger.gap2', '0.0081'),
                   ('uncharger.negativecharger.gap1', '0.0080'),
                   ('postivecharger.uncharger.gap5', '0.0080'),
                   ('aromatic.uncharger.gap0', '0.0080'), ('alphaticr.aromatic.gap2', '0.0080'),
                   ('alphaticr.uncharger.gap4', '0.0079'), ('uncharger.uncharger.gap2', '0.0079'),
                   ('alphaticr.postivecharger.gap3', '0.0078'),
                   ('uncharger.alphaticr.gap0', '0.0078'), ('aromatic.uncharger.gap5', '0.0077'),
                   ('uncharger.aromatic.gap1', '0.0076'),
                   ('postivecharger.alphaticr.gap1', '0.0076'),
                   ('uncharger.alphaticr.gap5', '0.0076'),
                   ('aromatic.negativecharger.gap0', '0.0076'),
                   ('negativecharger.alphaticr.gap5', '0.0076'),
                   ('aromatic.postivecharger.gap3', '0.0076'),
                   ('aromatic.aromatic.gap4', '0.0076'),
                   ('aromatic.negativecharger.gap5', '0.0075'),
                   ('alphaticr.negativecharger.gap4', '0.0075'),
                   ('aromatic.negativecharger.gap4', '0.0075'),
                   ('postivecharger.postivecharger.gap1', '0.0075'),
                   ('alphaticr.alphaticr.gap5', '0.0074'), ('alphaticr.uncharger.gap3', '0.0074'),
                   ('postivecharger.uncharger.gap0', '0.0074'),
                   ('negativecharger.uncharger.gap2', '0.0074'),
                   ('alphaticr.alphaticr.gap0', '0.0074'), ('alphaticr.alphaticr.gap1', '0.0074'),
                   ('uncharger.postivecharger.gap5', '0.0074'),
                   ('negativecharger.aromatic.gap2', '0.0073'),
                   ('aromatic.postivecharger.gap5', '0.0073'),
                   ('uncharger.postivecharger.gap1', '0.0073'),
                   ('postivecharger.uncharger.gap3', '0.0073'),
                   ('uncharger.aromatic.gap0', '0.0072'),
                   ('uncharger.postivecharger.gap2', '0.0072'),
                   ('uncharger.alphaticr.gap2', '0.0072'), ('uncharger.alphaticr.gap3', '0.0072'),
                   ('uncharger.alphaticr.gap1', '0.0072'), ('aromatic.aromatic.gap3', '0.0072'),
                   ('negativecharger.uncharger.gap3', '0.0072'),
                   ('postivecharger.alphaticr.gap4', '0.0071'),
                   ('uncharger.postivecharger.gap0', '0.0071'),
                   ('aromatic.negativecharger.gap3', '0.0071'),
                   ('alphaticr.postivecharger.gap4', '0.0071'),
                   ('alphaticr.aromatic.gap4', '0.0071'),
                   ('negativecharger.alphaticr.gap0', '0.0071'),
                   ('alphaticr.postivecharger.gap2', '0.0071'),
                   ('postivecharger.alphaticr.gap0', '0.0070'),
                   ('postivecharger.uncharger.gap2', '0.0070'),
                   ('alphaticr.aromatic.gap0', '0.0070'),
                   ('uncharger.postivecharger.gap4', '0.0070'),
                   ('aromatic.uncharger.gap1', '0.0070'),
                   ('negativecharger.negativecharger.gap3', '0.0069'),
                   ('uncharger.aromatic.gap4', '0.0069'), ('aromatic.uncharger.gap3', '0.0069'),
                   ('uncharger.negativecharger.gap3', '0.0069'),
                   ('uncharger.negativecharger.gap5', '0.0068'),
                   ('alphaticr.negativecharger.gap1', '0.0066'),
                   ('postivecharger.aromatic.gap5', '0.0065'),
                   ('negativecharger.aromatic.gap3', '0.0064'),
                   ('aromatic.negativecharger.gap1', '0.0063'),
                   ('alphaticr.uncharger.gap0', '0.0063'),
                   ('postivecharger.negativecharger.gap1', '0.0062'),
                   ('postivecharger.postivecharger.gap5', '0.0062'),
                   ('postivecharger.negativecharger.gap0', '0.0062'),
                   ('uncharger.postivecharger.gap3', '0.0062'),
                   ('postivecharger.negativecharger.gap2', '0.0062'),
                   ('alphaticr.uncharger.gap2', '0.0061'), ('aromatic.alphaticr.gap1', '0.0061'),
                   ('aromatic.postivecharger.gap2', '0.0061'),
                   ('aromatic.postivecharger.gap0', '0.0060'),
                   ('alphaticr.uncharger.gap5', '0.0060'), ('uncharger.uncharger.gap5', '0.0060'),
                   ('aromatic.alphaticr.gap5', '0.0060'), ('aromatic.aromatic.gap0', '0.0060'),
                   ('uncharger.negativecharger.gap0', '0.0060'),
                   ('negativecharger.alphaticr.gap1', '0.0059'),
                   ('alphaticr.negativecharger.gap0', '0.0059'),
                   ('uncharger.uncharger.gap4', '0.0059'),
                   ('negativecharger.negativecharger.gap5', '0.0059'),
                   ('negativecharger.negativecharger.gap0', '0.0059'),
                   ('postivecharger.postivecharger.gap2', '0.0059'),
                   ('postivecharger.negativecharger.gap4', '0.0058'),
                   ('alphaticr.postivecharger.gap1', '0.0058'),
                   ('postivecharger.aromatic.gap2', '0.0058'),
                   ('uncharger.negativecharger.gap4', '0.0058'),
                   ('alphaticr.alphaticr.gap2', '0.0058'), ('aromatic.aromatic.gap1', '0.0058'),
                   ('negativecharger.negativecharger.gap2', '0.0058'),
                   ('negativecharger.aromatic.gap4', '0.0058'),
                   ('uncharger.aromatic.gap2', '0.0058'), ('alphaticr.uncharger.gap1', '0.0058'),
                   ('negativecharger.uncharger.gap1', '0.0058'),
                   ('aromatic.aromatic.gap5', '0.0058'),
                   ('negativecharger.uncharger.gap4', '0.0058'),
                   ('alphaticr.aromatic.gap1', '0.0058'), ('uncharger.uncharger.gap0', '0.0058'),
                   ('uncharger.alphaticr.gap4', '0.0057'),
                   ('negativecharger.uncharger.gap0', '0.0057'),
                   ('negativecharger.postivecharger.gap3', '0.0057'),
                   ('alphaticr.postivecharger.gap5', '0.0057'),
                   ('negativecharger.uncharger.gap5', '0.0057'),
                   ('negativecharger.postivecharger.gap0', '0.0057'),
                   ('negativecharger.alphaticr.gap2', '0.0057'),
                   ('postivecharger.postivecharger.gap0', '0.0057'),
                   ('negativecharger.aromatic.gap0', '0.0057'),
                   ('negativecharger.aromatic.gap5', '0.0057'),
                   ('uncharger.uncharger.gap1', '0.0057'),
                   ('negativecharger.negativecharger.gap4', '0.0057'),
                   ('postivecharger.aromatic.gap0', '0.0057'),
                   ('negativecharger.postivecharger.gap1', '0.0057'),
                   ('alphaticr.aromatic.gap5', '0.0057'),
                   ('postivecharger.aromatic.gap1', '0.0057'),
                   ('uncharger.negativecharger.gap2', '0.0057'),
                   ('alphaticr.alphaticr.gap3', '0.0056'), ('aromatic.aromatic.gap2', '0.0056'),
                   ('aromatic.negativecharger.gap2', '0.0056'),
                   ('negativecharger.aromatic.gap1', '0.0056'),
                   ('postivecharger.postivecharger.gap3', '0.0056'),
                   ('postivecharger.alphaticr.gap2', '0.0056'),
                   ('negativecharger.postivecharger.gap4', '0.0056'),
                   ('alphaticr.aromatic.gap3', '0.0056'),
                   ('aromatic.postivecharger.gap1', '0.0055'),
                   ('postivecharger.alphaticr.gap5', '0.0055'),
                   ('negativecharger.negativecharger.gap1', '0.0055'),
                   ('negativecharger.postivecharger.gap2', '0.0055'),
                   ('postivecharger.uncharger.gap1', '0.0055'),
                   ('negativecharger.postivecharger.gap5', '0.0055'),
                   ('alphaticr.negativecharger.gap2', '0.0055'),
                   ('aromatic.alphaticr.gap3', '0.0054'),
                   ('alphaticr.postivecharger.gap0', '0.0053'),
                   ('postivecharger.aromatic.gap4', '0.0052'),
                   ('aromatic.alphaticr.gap2', '0.0050'), ('aromatic.alphaticr.gap0', '0.0042'),
                   ('postivecharger.negativecharger.gap3', '0.0041')]

    def metrics_Line1(metrics_file):
        df = pd.read_csv(metrics_file)
        x_axis = df['length']
        df.pop('length')
        line = Line("metrics line")
        max_v = max(df.max())
        min_v = min(df.min())
        for metric_name in df:
            line.add(
                metric_name,
                x_axis,
                df[metric_name],
                yaxis_max=float("{0:.2f}".format(max_v + max_v * 0.1)),
                yaxis_min=float("{0:.2f}".format(min_v - min_v * 0.1)),
                xaxis_type='value'
            )
        # print(max_v+max_v*0.1)
        # print(min_v - min_v*0.1)
        # clean_upload_file(os.path.dirname(metrics_file))  # delete Results folders
        return line

    line = metrics_Line1(
        os.path.dirname(os.path.abspath(__file__)) + os.sep + 'static' + os.sep + '150.arff.metrics.csv')

    return  render_template('result.html', title='result', accuarcy='0.7945', precision='0.8005', recall='0.8000', f1='0.7992',
                           roc_area='0.7963',
                           TP='50', FN='38', FP='13', TN='9',
                           zip_name='',
                           myechart=line.render_embed(),
                           host=REMOTE_HOST,
                           script_list=line.get_js_dependencies(),
                           features_kv=features_kv,
                           unique='150_example.zip',
                           download = True)


if __name__ == '__main__':
    app.run(debug=False,host='0.0.0.0',port='10503',threaded=True)


