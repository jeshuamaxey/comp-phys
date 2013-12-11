'use strict';

var fs = require('fs');

//deal with that pesky .DS_Store file in OSX
//fs.unlink('./videos/.DS_Store', startUpdating(), function(err){});

function main() {
	//get dir listing as an array
	walk('./data', outputToFile); //fs.readdirSync('./data');
}

//taken from incredible stackoverflow answer
//http://stackoverflow.com/questions/5827612/node-js-fs-readdir-recursive-directory-search
function walk(dir, done) {
  var results = [];
  fs.readdir(dir, function(err, list) {
    if (err) return done(err);
    var i = 0;
    (function next() {
      var file = list[i++];
      if (!file) return done(null, results);
      file = dir + '/' + file;
      fs.stat(file, function(err, stat) {
        if (stat && stat.isDirectory()) {
          walk(file, function(err, res) {
            results = results.concat(res);
            next();
          });
        } else {
          results.push(file);
          next();
        }
      });
    })();
  });
};

function outputToFile(err, originalListing) {

	var str = '';
	originalListing.forEach(function(file) {
		str += file + ',';
	});

	//write str to file
	var outputFile = './data/data-files.csv';

	fs.writeFile(outputFile, str, function(err) {
		if(err) {
			console.log(err);
		} else {
			console.log("Files listed in " + outputFile + "\n\n");
		}
	});
	return false;
}

//GO!
main();
