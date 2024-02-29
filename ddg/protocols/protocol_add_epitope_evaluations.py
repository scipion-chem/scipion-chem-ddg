# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

from pwem.protocols import EMProtocol
from pyworkflow.protocol import params

from pwchem.objects import SetOfSequenceROIs

from .. import Plugin as ddgPlugin
from ..constants import EVAL_PARAM_MAP

class ProtDDGEvaluations(EMProtocol):
  """Run evaluations on a set of epitopes (SetOfSequenceROIs)"""
  _label = 'ddg epitope evaluations'

  _evaluatorOptions = ['Vaxijen2', 'Vaxijen3', 'AllerTop2', 'AllergenFP1']

  _vaxiTargets = ['bacteria', 'virus', 'tumor', 'parasite', 'fungal']

  _softParams = {'Vaxijen2': ['vaxi2Target'],
                 'Vaxijen3': [],
                 'AllerTop': [],
                 'AllergenFP': [],
                 }

  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)

  def _defineParams(self, form):
    form.addSection(label='Input')
    iGroup = form.addGroup('Input')
    iGroup.addParam('inputROIs', params.PointerParam, pointerClass="SetOfSequenceROIs", label='Input epitopes: ',
                    help="Input set of epitope sequences as SetOfSequenceROIs")

    form.addSection(label='Add evaluations')
    aGroup = form.addGroup('Define evaluator')
    aGroup.addParam('chooseEvaluator', params.EnumParam, choices=self._evaluatorOptions,
                    label='Choose evaluator: ', default=0,
                    help='Epitope evaluation software to use.')
    aGroup.addParam('evaluatorName', params.StringParam, label='Evaluator name: ',
                    default='', expertLevel=params.LEVEL_ADVANCED,
                    help='Set the name for the defined evaluator.')

    aGroup.addParam('vaxi2Target', params.EnumParam, choices=self._vaxiTargets, default=0,
                    label='Vaxijen2 target: ', condition='chooseEvaluator==0',
                    help='Target type for the Vaxijen2 epitopen evaluation')

    aGroup.addParam('addEval', params.LabelParam, label='Add defined evaluator: ',
                    help='Add defined evaluator to perform the epitope prediction')

    sGroup = form.addGroup('Evaluators summary')
    sGroup.addParam('inEvals', params.TextParam, width=70, default='',
                    label='Evaluators summary: ',
                    help='Summary of the epitope evaluations that will be performed')

    form.addParallelSection(threads=4, mpi=1)


  def _insertAllSteps(self):
    self._insertFunctionStep(self.evaluationStep)

  def evaluationStep(self):
    nt = self.numberOfThreads.get()
    sDics = self.getWebEvaluatorDics()
    sequences = self.getInputSequences()

    epiDic = ddgPlugin.performEvaluations(sequences, sDics, nt, ddgPlugin.getBrowserData())
    print(epiDic)

    outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
    for i, roi in enumerate(self.inputROIs.get()):
      for (evalKey, softName), scores in epiDic.items():
        setattr(roi, evalKey, params.Float(scores[i]))
      outROIs.append(roi)

    if len(outROIs) > 0:
      self._defineOutputs(outputROIs=outROIs)


  ##################### UTILS #####################
  def getInputSequences(self):
    seqs = {}
    for roi in self.inputROIs.get():
      seqs[roi.getROIId()] = roi.getROISequence()
    return seqs

  def buildElementDic(self):
    sName, soft = self.evaluatorName.get(), self.getEnumText('chooseEvaluator')
    if not sName.strip():
      sName = self.getDefSName(soft)

    sDic = {sName: {'software': soft}}
    for paramName in self._softParams[soft]:
      sDic[sName].update({paramName: self.getParamValue(paramName)})
    return sDic

  def parseElementsDic(self):
    ''' Parse the selector dictionaries included in the input list
    :return: dic, {selName: {software: softName, paramName: paramValue}} with the chosen Scipion parameters
    '''
    sDic = {}
    for line in self.inEvals.get().split('\n'):
      if line.strip():
        sd = f'{{{line.split(") ")[1]}}}'
        sDic.update(eval(sd))
    return sDic

  def getDefSName(self, soft):
    sDic, i = self.parseElementsDic(), 1
    sName = f'{soft}-{i}'
    while sName in sDic:
      i += 1
      sName = f'{soft}-{i}'
    return sName

  def getParamValue(self, paramName):
    if isinstance(self.getParam(paramName), params.EnumParam):
      value = self.getEnumText(paramName)
    else:
      value = getattr(self, paramName).get()
    return value

  def getWebEvaluatorDics(self):
    ''' Returns the selector dictionary with the parameter names expected by the web server
    :return: dic, {selName: {software: softName, paramName: paramValue}} with the webserver chosen parameters
    '''
    wsDic = {}
    sDic = self.parseElementsDic()
    for sName, curSDic in sDic.items():
      wsDic[sName] = {}
      for paramName, paramValue in curSDic.items():
        if paramName in EVAL_PARAM_MAP:
          paramName = EVAL_PARAM_MAP[paramName]
        wsDic[sName][paramName] = paramValue
    return wsDic