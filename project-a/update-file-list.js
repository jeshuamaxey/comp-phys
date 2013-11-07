'use strict';

var fs = require('fs');

//deal with that pesky .DS_Store file in OSX
//fs.unlink('./videos/.DS_Store', startUpdating(), function(err){});

function main() {

	//get dir listing as an array
	var originalListing = fs.readdirSync('./data');

	var str = '';

	originalListing.forEach(function(file) {
		str += 'data/' + file + ',';
	});

	//write str to file
	var outputFile = './data/data-files.csv';

	fs.writeFile(outputFile, str, function(err) {
		if(err) {
			console.log(err);
		} else {
			console.log("Files listed in " + outputFile);
		}
	});
	return false;
}

main();
