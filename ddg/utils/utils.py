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

import time, requests
from Bio import SeqIO

from ..constants import EVAL_PARAM_MAP

def parseInputProteins(faFile):
  '''Uses BioPython to parse a fasta file and return it as dictionary
  :param faFile: input fasta filename
  :return: {seqName1: seqStr1, ...}
  '''
  faDic = {}
  with open(faFile) as handle:
    for values in SeqIO.FastaIO.SimpleFastaParser(handle):
        faDic[values[0]] = values[1]
  return faDic


def reportPoolStatus(poolDic):
  '''Check the status of the AsynPool objects stored as values of the dictionary and reports when they finish
  '''
  ready = []
  while len(ready) < len(poolDic):
    time.sleep(5)
    for evalSoft, po in poolDic.items():
      if po.ready() and evalSoft not in ready:
        ready.append(evalSoft)
        print(f'{evalSoft} execution finished ({len(ready)} / {len(poolDic)})')


def divide_chunks(iter, chunkSize):
  '''Divides an iterable into chunks of size chunkSize'''
  chunks = []
  for i in range(0, len(iter), chunkSize):
    chunks.append(iter[i:i + chunkSize])
  return chunks

def buildSeqFasta(seqLists):
  '''From a list of sequence chunks, build a list of those sequences fasta strings'''
  seqStrs = []
  for seqList in seqLists:
    fastaList = [f'>seq{i+1}\n{seq}\n' for i, seq in enumerate(seqList)]
    seqStrs.append(''.join(fastaList))
  return seqStrs

def getFastaStrs(seqDic, maxChunk=None):
  '''Build a list of fasta strings from a list of sequences in chunks of maxChunk size'''
  maxChunk = len(seqDic) if not maxChunk else maxChunk
  seqList = list(seqDic.values())
  seqLists = divide_chunks(seqList, maxChunk)
  return buildSeqFasta(seqLists)

def getFastaFiles(seqDic, evalSoft, maxChunk=None):
  '''Write a series of fasta files with maxChunk number of sequences from a set of sequences'''
  fastaStrs = getFastaStrs(seqDic, maxChunk)
  faFiles = []
  for i, fStr in enumerate(fastaStrs):
    faFiles.append(f'/tmp/{evalSoft}_input_{i}.fa')
    with open(faFiles[-1], 'w') as f:
      f.write(fStr)
  return faFiles

def setData(driver, paramDic):
  '''Sets the additional data parameters in the web of the software evaluation
  driver: selenium driver, with url set in the software web
  paramDic: dic, containing the parameter names (as dic keys) and corresponding values (dic values)
  '''
  from selenium.webdriver.common.by import By
  for dk, dv in paramDic.items():
    dataElements = driver.find_elements(By.NAME, dk)
    for dEl in dataElements:
      if dEl.get_attribute('value') == dv:
        dEl.click()
  return driver


def getDriver(browserData):
  from selenium import webdriver
  from selenium.webdriver.chrome.options import Options as ChromeOptions
  from selenium.webdriver.firefox.options import Options as FireOptions
  '''Return a selenium WebDriver object, depending on the selected browser
  - browserData: dic, contains the information about the browser to be used
    - name: str, the name of the browser to use (either "Chrome" for Google-Chrome or Firefox)
    - path: str, path for the browser executable in case of non default
  '''
  if not 'name' in browserData or browserData['name'] != 'Firefox':
    options = ChromeOptions()
    driverObj = webdriver.Chrome
    browserPath = '/usr/bin/google-chrome' if (not 'path' in browserData or not browserData['path'])\
      else browserData['path']
  else:
    options = FireOptions()
    driverObj = webdriver.Firefox
    browserPath = '/usr/bin/firefox' if (not 'path' in browserData or not browserData['path']) \
      else browserData['path']

  options._binary_location = browserPath
  options.add_argument('--headless')
  driver = driverObj(options=options)
  return driver


def performRequest(seqKeys, driver, softData):
  from selenium.webdriver.common.by import By
  '''Performs a request in a evaluation software using selenium to emulate the browser.
  - seqData: dic, contains the keys and values of the web elements to write, including the
             sequence in the format expected by the web (fasta file, fasta string or sequence string)
  - driver: selenium driver to use for the request
  - softData: dic, containing all the characteristics and info for the specific sofware web. Among others (key: value):
    - url: str, evaluation software url
    - params: dict, additional data arguments to fill in the web form
    - submitCSS: str, css selector to identify the submit button (e.g: "input[name='Submit']")
  - seqKeys: dic, if not None, specifies the web html name key and value to write the sequence name. e.g: {seqName: seq1}
  '''
  driver.get(softData['url'])

  for xKeyName, xKeyVal in seqKeys.items():
    extraElem = driver.find_element(By.NAME, xKeyName)
    extraElem.send_keys(xKeyVal)

  driver = setData(driver, softData['params'])
  driver.find_elements(By.CSS_SELECTOR, softData['submitCSS'])[0].click()
  return driver


def getSeqData(seqDic, softData):
  '''Returns a list containing the chunks of sequences as expected from the web to use.
  It can be either: a list with one fasta file, a list with one fasta string or a list with sequences strings
  - seqDic: dic, sequences {seqId: seqString}
  - softData: dic, containing all the characteristics and info for the specific sofware web. Among others (key: value):
    - multi: whether the web admits multiple sequences at one time
    - seqFormat: whether to return a fasta file ("fastaFile") or the fasta string ("fastaString")
    - softName: software name for the fasta file to be named
  '''
  if softData['multi']:
    if softData['seqFormat'] == 'fastaFile':
      seqData = getFastaFiles(seqDic, softData['softName'])
    else:
      seqData = getFastaStrs(seqDic)
  else:
    seqData = list(seqDic.values())
  return seqData


def updateBatchDic(outDic, batchDic):
  '''Updates(appends) the lists inside the outDic values with the ones in the batchDic
  '''
  for key, values in batchDic.items():
    if key in outDic:
      outDic[key] += values
    else:
      outDic[key] = values
  return outDic


def seleniumRequest(seqDic, softData, browserData, parseFunction, seqNameKey=None):
  '''Perform a series of Selenium requests an operations to emulate the evaluation of a set of sequences by a software
  web server.
  - seqDic: dic, sequences {seqId: seqString}
  - softData: dic, contains the information necessary to build the software web request
  - browserData: dic, contains the information necessary to build the Selenium driver
  - parseFunction: func, parses the driver data once the request is performed and returns a dic {'Score' [sc1, ...]}
  - seqNameKey: str, if not None, include the sequence name as a web element value to write in this key
  '''
  # url, data, softName, seqFormat='fastaString', seqName='sequence', multi=True
  driver = getDriver(browserData)
  seqData = getSeqData(seqDic, softData)

  # Performing one request for each chunk of admitted data (just once if fasta admitted)
  outDic = {}
  for i, seq in enumerate(seqData):
    curSeqKeys = {softData['seqName']: seq}
    if seqNameKey:
      curSeqKeys.update({seqNameKey: f'seq{i + 1}'})

    driver = performRequest(curSeqKeys, driver, softData)
    # Parse the driver with the corresponding function for each software
    batchDic = parseFunction(driver)
    outDic = updateBatchDic(outDic, batchDic)
  return outDic


def innerSplit(text, preText, endText):
  results, splitted = [], text.split(preText)[1:]
  for text in splitted:
    results.append(text.strip().split(endText)[0].strip())
  return results

########## REQUESTS ##########

def makeRequest(url, action='post', data={}, headers={}):
  if action == 'post':
    response = requests.post(url, data=data, headers=headers)
  else:
    response = requests.get(url, data=data, headers=headers)

  if response.status_code == 200:
    pass
  else:
    print(f"There was an error in request to {url}: {response.status_code}")
  return response


########### SELENIUM CALLS ################

def callVaxijen3(sequences, browserData={}, data={}):
  softData = {'url': "https://www.ddg-pharmfac.net/vaxijen3/",
              'multi': True, 'seqFormat': 'fastaFile', 'softName': 'Vaxijen3',
              'seqName': 'uploaded_file', 'params': data, 'submitCSS': "input[name='submit']"}

  outDic = seleniumRequest(sequences, softData, browserData, parseVaxijen3)
  return outDic


def callVaxijen2(sequences, browserData={}, data={}):
  data = {"Target": 'Bacteria'} if not data else data

  softData = {'url': "https://www.ddg-pharmfac.net/vaxijen/VaxiJen/VaxiJen.html",
              'multi': True, 'seqFormat': 'fastaFile', 'softName': 'Vaxijen2',
              'seqName': 'uploaded_file', 'params': data, 'submitCSS': "input[name='submit']"}

  outDic = seleniumRequest(sequences, softData, browserData, parseVaxijen2)

  return outDic


def callAllerTop2(sequences, browserData={}, data={}):
  softData = {'url': "https://www.ddg-pharmfac.net/AllerTOP/",
              'multi': False,
              'seqName': 'sequence', 'params': data, 'submitCSS': "input[name='Submit']"}

  outDic = seleniumRequest(sequences, softData, browserData, parseAllerDDG)

  return outDic


def callAllergenFP1(sequences, browserData={}, data={}):
  softData = {'url': "https://ddg-pharmfac.net/AllergenFP/",
              'multi': False,
              'seqName': 'sequence', 'params': data, 'submitCSS': "input[name='Submit']"}

  outDic = seleniumRequest(sequences, softData, browserData, parseAllerDDG)
  return outDic

############## PARSING ##############

def parseVaxijen3(driver):
  from selenium.webdriver.common.by import By
  data = driver.find_elements(By.CSS_SELECTOR, "table[class='boilerplate']")
  while not data:
    time.sleep(5)
    data = driver.find_elements(By.CSS_SELECTOR, "table[class='boilerplate']")
  resultText = data[0].text

  results = innerSplit(resultText, 'is predicted to be', 'with')
  probs = innerSplit(resultText, 'with probability', '\n')

  resDic = {'Score': []}
  for i in range(len(results)):
    posNeg = 1 if results[i] == 'Probable ANTIGEN' else -1
    resDic['Score'].append(float(probs[i].replace('%', '')) * posNeg * 0.01)

  return resDic


def parseVaxijen2(driver):
  from selenium.webdriver.common.by import By
  data = driver.find_elements(By.CSS_SELECTOR, "table[border='0']")
  while not data:
    time.sleep(5)
    data = driver.find_elements(By.CSS_SELECTOR, "table[border='0']")
  resultText = data[0].text

  results = innerSplit(resultText, '(', ')')
  probs = innerSplit(resultText, '=', '(')

  resDic = {'Score': []}
  for i in range(len(results)):
    posNeg = 1 if results[i] == 'Probable ANTIGEN' else -1
    resDic['Score'].append(float(probs[i].replace('%', '')) * posNeg)

  return resDic


def parseAllerDDG(driver):
  from selenium.webdriver.common.by import By
  data = driver.find_elements(By.CSS_SELECTOR, "table[border='0']")
  while not data:
    time.sleep(5)
    data = driver.find_elements(By.CSS_SELECTOR, "table[border='0']")
  resultText = data[0].text

  results = innerSplit(resultText, 'Your sequence is:\n', '\n')

  resDic = {'Score': []}
  for i in range(len(results)):
    resDic['Score'].append(0 if results[i] == 'PROBABLE NON-ALLERGEN' else 1)
  return resDic

def mapEvalParamNames(sDic):
  wsDic = {}
  for sName, curSDic in sDic.items():
    wsDic[sName] = {}
    for paramName, paramValue in curSDic.items():
      if paramName in EVAL_PARAM_MAP:
        paramName = EVAL_PARAM_MAP[paramName]
      wsDic[sName][paramName] = paramValue
  return wsDic
