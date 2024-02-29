# **************************************************************************
# *
# * Authors:	Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This package contains protocols for creating and using IIITD Raghava software
"""

import multiprocessing

from scipion.install.funcs import InstallHelper

from pwchem import Plugin as pwchemPlugin

from .utils import *
from .constants import *

# Pluging variables
_logo = 'ddg_logo.png'

class Plugin(pwchemPlugin):
	"""
	"""

	@classmethod
	def _defineVariables(cls):
		cls._defineVar(DDG_DIC['browser'], 'Chrome')
		cls._defineVar(DDG_DIC['browserPath'], '/usr/bin/google-chrome')

	@classmethod
	def defineBinaries(cls, env, default=True):
		"""This function defines the binaries for each package."""
		installer = InstallHelper(DDG_DIC['name'], packageHome=cls.getVar(DDG_DIC['home']),
															packageVersion=DDG_DIC['version'])

		# Installing package
		installer.addPackage(env, ['git'], default=default)

	# ---------------------------------- Protocol functions-----------------------
	@classmethod
	def performEvaluations(cls, sequences, evalDics, jobs=1, browserData={}):
		'''Generalize caller to the evaluation functions.
    - sequences: list of sequences
    - evalDics: dictionary as {(evalKey, softwareName): {parameterName: parameterValue}}
    - jobs: int, number of jobs for parallelization
    Returns a dictionary of the form: {(evalKey, softwareName): {'Score': [scoreValues], otherParam: [otherValues]}}
    '''
		funcDic = {
			'Vaxijen2': callVaxijen2, 'Vaxijen3': callVaxijen3, 'AllerTop2': callAllerTop2, 'AllergenFP1': callAllergenFP1
		}

		# Create a pool of worker processes
		nJobs = len(evalDics) if len(evalDics) < jobs else jobs
		pool = multiprocessing.Pool(processes=nJobs)

		resultsDic = {}
		for evalKey, evalDic in evalDics.items():
			softName = evalDic['software']
			smallEvalDic = evalDic.copy()
			del smallEvalDic['software']
			if softName in funcDic:
				resultsDic[(evalKey, softName)] = pool.apply_async(funcDic[softName],
																													 args=(sequences, smallEvalDic, browserData))

		reportPoolStatus(resultsDic)

		pool.close()
		pool.join()

		epiDics = {}
		for (evalKey, softName), res in resultsDic.items():
			epiDics[(evalKey, softName)] = res.get()['Score']

		return epiDics

	# ---------------------------------- Utils functions-----------------------
	@classmethod
	def getBrowserData(cls):
		return {'name': cls.getVar(DDG_DIC['browser']), 'path': cls.getVar(DDG_DIC['browserPath'])}