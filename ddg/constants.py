# coding: latin-1
# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

# Common constants
DEFAULT_VERSION = '3.0'

# Package dictionaries
DDG_DIC = {'name': 'DDG',    'version': '3.0',
           'home': 'DDG_HOME', 'activation': 'DDG_ACTIVATION_CMD',
           'browser': 'DDG_BROWSER', 'browserPath': 'DDG_BROWSER_PATH'}

EVAL_PARAM_MAP = {'ToxinPred': {'method': {'SVM (Swiss-Prot)': 1, 'SVM (Swiss-Prot) + Motif': 2, 'SVM (TrEMBL)': 3}}}

EVALSUM = '''1) "Vaxijen2-1": {'software': 'Vaxijen2', 'vaxi2Target': 'bacteria'}
2) "Vaxijen3-1": {'software': 'Vaxijen3'}
3) "AllerTop2-1": {'software': 'AllerTop2'}
4) "AllergenFP1-1": {'software': 'AllergenFP1'}'''