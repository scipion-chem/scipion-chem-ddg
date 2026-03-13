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
from ..utils import mapEvalParamNames

class ProtDDGEvaluations(EMProtocol):
  """Run evaluations on a set of epitopes (SetOfSequenceROIs)

    AI Generated:

      ProtDDGEvaluations - User Manual

      Overview
      --------
      The AddEpitopeEvaluations protocol enables the annotation of protein mutations
      with epitope-related information. It integrates curated or predicted epitope data
      into mutation-level analyses, supporting studies on antigen–antibody interactions,
      immune escape, or immunogenicity assessment. This protocol is particularly useful
      for workflows in structure-guided immunology, vaccine design, or antibody engineering
      within the Scipion-Chem environment.

      Input Requirements
      ------------------
      1. **Mutation Dataset**:
         - Provide a `SetOfSequenceROIs` object containing mutations or sequence regions.
         - Each mutation should have residue-level detail and reference sequence context.

      2. **Evaluator Selection**:
         - Choose one or more evaluation tools from:
           - Vaxijen2
           - Vaxijen3
           - AllerTop2
           - AllergenFP1
         - Configure evaluator-specific parameters such as Vaxijen2 target type (e.g., bacteria, virus, tumor, parasite, fungal).

      3. **Evaluator Metadata** (optional):
         - Define an evaluator name for identification in outputs.
         - Summarize configured evaluators for review before execution.

      Workflow
      --------
      1. **Parameter Definition**:
         - Specify the evaluation software and associated settings.
         - Add multiple evaluators sequentially if needed.

      2. **Input Parsing**:
         - Extract mutation sequences from the input `SetOfSequenceROIs`.
         - Map mutation positions to known or predicted epitope regions.

      3. **Epitope Evaluation**:
         - Execute the selected evaluation tools.
         - Assign scores or flags to each mutation based on epitope presence or likelihood.
         - Thresholds and scoring options can be applied to distinguish core versus peripheral epitopes.

      4. **Output Generation**:
         - Annotate the input dataset with new columns reflecting epitope evaluations.
         - Save results as an enriched `SetOfSequenceROIs` object.
         - Outputs can be visualized, filtered, or exported for downstream analysis.

      Outputs
      -------
      - **Annotated Mutations**: `SetOfSequenceROIs` with epitope evaluation scores or flags for each mutation.
      - **Evaluator Summary**: Text summary of the evaluators applied and their settings.

      Validation & Warnings
      ---------------------
      - Ensure the input dataset contains residue-level detail for proper mapping.
      - Check that the evaluator parameters match the chosen software requirements.
      - Sequential evaluator addition is recommended to prevent conflicts.
      - The protocol does not modify input sequences; it only adds annotation columns.

      Practical Recommendations
      -------------------------
      - Start with high-confidence mutation or sequence datasets.
      - Use multiple evaluators for complementary epitope information.
      - Verify evaluator targets (e.g., bacteria vs. virus) are appropriate for your study.
      - Inspect enriched datasets visually or statistically to validate predicted epitope impact.

      Final Perspective
      -----------------
      AddEpitopeEvaluations provides a reproducible method to enrich mutation-level
      datasets with immunological context. It integrates epitope mapping information
      from experimental or computational sources, supports flexible evaluator
      configuration, and produces outputs suitable for visualization, prioritization,
      and downstream immunological analyses in Scipion-Chem.
  """
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

  def _defineEvalParams(self, aGroup, allCond=True):
    '''Define the evaluation options and the parameters for each of them.
    allCond: condition to apply for all the parameters

    WARNING: This function is used by a scipion-chem metaprotocol to use and define this parameters by its own,
    modify with care
    '''
    aGroup.addParam('chooseDDGEvaluator', params.EnumParam, choices=self._evaluatorOptions,
                    label='Choose evaluator: ', default=0, condition=f'{allCond}',
                    help='Epitope evaluation software to use.')

    aGroup.addParam('vaxi2Target', params.EnumParam, choices=self._vaxiTargets, default=0,
                    label='Vaxijen2 target: ', condition=f'{allCond} and chooseDDGEvaluator==0',
                    help='Target type for the Vaxijen2 epitopen evaluation')
    return aGroup

  def _defineParams(self, form):
    form.addSection(label='Input')
    iGroup = form.addGroup('Input')
    iGroup.addParam('inputROIs', params.PointerParam, pointerClass="SetOfSequenceROIs", label='Input epitopes: ',
                    help="Input set of epitope sequences as SetOfSequenceROIs")

    form.addSection(label='Add evaluations')
    aGroup = form.addGroup('Define evaluator')
    aGroup = self._defineEvalParams(aGroup)
    aGroup.addParam('evaluatorDDGName', params.StringParam, label='Evaluator name: ',
                    default='', expertLevel=params.LEVEL_ADVANCED,
                    help='Set the name for the defined evaluator.')
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
    sName, soft = self.evaluatorName.get(), self.getEnumText('chooseDDGEvaluator')
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
    sDic = self.parseElementsDic()
    wsDic = mapEvalParamNames(sDic)
    return wsDic