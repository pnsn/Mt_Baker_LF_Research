// Mt. Baker Analyst Review of Well Located Events
// Script Generated from ChatGPT 3.5 Prompt by N. Stevens
// Script adapted & tested by N. Stevens for functionality before distribution (WIP)
function createClassificationForm() {
  // Open the Google Form
  var form = FormApp.create('Mt. Baker Event Classification Review')
    .setDescription('This form contains the EVIDs of "well-located" seismic events within 20km of the summit of Mt. Baker. Please help by reviewing them in Jiggle and providing your assessment of what event type you think they are in this form.\n\nNotes: \n - Event type codes are described at the top of each page.\n - An ASCII table with the EVIDs that can be loaded into Jiggle is linked at the top of each page for your conveineice.\n - The EVIDs are purpously shuffled in this form to shake-up your perception, please try to work through them in the order listed, rather than sequential order.\n - You can edit your responses until the form is closed. Think of the "submit" button as a "save" button!\n\nAI Disclosure: Form programmatically generated using Google App Script and code examples from ChatGPT. Tested, modified, and validated by N. Stevens.');

  // Get active spreadsheet
  var sheet = SpreadsheetApp.getActiveSpreadsheet();

  // Access EVID sheet entries
  var evidSheet = sheet.getSheets()[0];
  var evidValues = evidSheet.getRange(1, 1, evidSheet.getLastRow()).getValues();

  // Shuffle EVID values
  evidValues = shuffleArray(evidValues);

  // Access Class sheet entries
  var classSheet = sheet.getSheets()[1];
  var classValues = classSheet.getRange(1,1,classSheet.getLastRow()).getValues();

  // Covert classification to 1D array
  var classOptions = classValues.map(function(row) {return row[0];
  });

  //set number of entries per page
  var evidsPerPage = 10;

  // Loop over EVIDs in batches of 10 to create pages
  for (var _e = 0; _e < evidValues.length; _e += evidsPerPage) {
    form.addPageBreakItem()
      .setTitle('Page ' + ((_e / evidsPerPage) + 1) + ' of ' + Math.ceil(evidValues.length/evidsPerPage))
      .setHelpText('Event Classes:\n\neq = earthquake\nlf eq = low frequency earthquake\nlf su = low frequency surface event\nsu = surface event\npx = probable blast\nuk = unknown\nother = exotic events (leave a comment please!)\n\nEVID ASCII Table Here:\nhttps://drive.google.com/file/d/1b-Rmsehd4mMJeP7w-ioCn6uP8hGrfarJ/view?usp=share_link');

    var gridItem = form.addGridItem()
      .setTitle('Classify the following catalog events:')
      .setColumns(classOptions);

      // Add up to 10 EVIDs (rows) in the grid
      var evidRows = [];
      for (var _f = _e; _f < _e + evidsPerPage && _f < evidValues.length; _f++) {
        evidRows.push(evidValues[_f][0]);
      }
    gridItem.setRows(evidRows);

    // add optional comment field for this page
    form.addParagraphTextItem()
      .setTitle('Optional: Add any comments on events for this page')
      .setHelpText('Comment Format: start with `EVID#, ` and end with a semicolon.\n\ne.g.,\n60135818, event in coda of earlier event (already in event remarks);\n60135817, Nate made this example evid comment;');
      
  }

  // // Loope over EVIDs and add a question to the form
  // for (var _e = 0; _e < evidValues.length; _e++) {
  //   var evidNumber = evidValues[_e][0];
  //   form.addMultipleChoiceItem()
  //   .setTitle('Event ID: ' + evidNumber)
  //   .setChoiceValues(classOptions);
  // }
  Logger.log('Form created: ' + form.getEditUrl());
}

function shuffleArray(array) {
  for (let _e = array.length - 1; _e > 0; _e--) {
    const _f = Math.floor(Math.random() * (_e + 1));
    [array[_e], array[_f]] = [array[_f], array[_e]];
  }
  return array;
}