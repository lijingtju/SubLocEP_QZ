from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField, SubmitField, BooleanField, TextAreaField, RadioField
from flask_wtf.file import FileField, FileAllowed
from wtforms.validators import DataRequired

class DataForm(FlaskForm):
    start_index = StringField('start_index',)
    end_index =  StringField('end_index')
    #length =  StringField('length', validators=[DataRequired()])
    step_length = StringField('length')
    EP = RadioField('EP', choices=[('EP1', 'EP1'), ('EP2', 'EP2')],
                                  validators=[DataRequired()])
    dataset = FileField('Upload Dataset', validators=[FileAllowed(['arff', 'libsvm','csv'])])
    submit = SubmitField('Start EP3')
