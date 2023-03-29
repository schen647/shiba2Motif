// start a websocket server on port 8080
const WebSocketServer = require("ws").Server;
const wss = new WebSocketServer({ port: 8089 });
const fs = require('fs');
const http = require('http');
const formidable = require('formidable');

const { spawn } = require("child_process");
const refList = [];
// get all files under predefinedref folder
fs.readdirSync('./predefinedref').forEach(file => {
  refList.push(file);
});


wss.on("connection", (ws) => {
  // create input$ip and ref$ip in tmp for each client
  const ip = ws._socket.remoteAddress;

  if (!fs.existsSync(`/tmp/input${ip}/`)) {
    fs.mkdirSync(`/tmp/input${ip}/`);
  }
  if (!fs.existsSync(`/tmp/ref${ip}`)) {
    fs.mkdirSync(`/tmp/ref${ip}`);
  }
  ws.reffa = `/tmp/ref${ip}/ref.fa`;
  ws.on("message", (message) => {
    // parse the message as json
    const data = JSON.parse(message);
    const action = data.action;
    const parameters = data.parameters;
    const ip = ws._socket.remoteAddress;
    switch (action) {

      case "run":
        console.log('launching worker1')
        // if the ws already has a python script runing, return it
        if (ws.py)
          return;


        const tmpDir = `/tmp/meme${ip}`;
        const out = `/tmp/out${ip}`;
        const inputDir = `/tmp/input${ip}/`;
        //const ref = ws.reffa;
        console.log('launching worker')
        console.log("./worker.py",
        '--tmp', "'"+tmpDir+"'",
        '--out', "'"+out+"'",
        '--input_five', "'"+inputDir + 'five.txt'+"'",
        '--input_mxe', "'"+inputDir + 'mxe.txt'+"'",
        '--input_ri', "'"+inputDir + 'ri.txt'+"'",
        '--input_se', "'"+inputDir + 'se.txt'+"'",
        '--input_three', "'"+inputDir + 'three.txt'+"'",
        '--ref', "'"+ws.reffa+"'")
        ws.py = spawn("/usr/bin/python3", ["./worker.py",
          '--tmp', tmpDir,
          '--out', out,
          '--input_five', inputDir + 'five.txt',
          '--input_mxe', inputDir + 'mxe.txt',
          '--input_ri', inputDir + 'ri.txt',
          '--input_se', inputDir + 'se.txt',
          '--input_three', inputDir + 'three.txt',
          '--ref', ws.reffa,
          '--no_detach'
        ]);
        // on return, send the data to the client
        ws.py.stdout.on("data", (data) => {
          ws.send(JSON.stringify({ action: 'data', parameters: { 'output': data } }));
          ws.py = null;
        });
        break;

      case "ident":
        console.log("backend ack");
        ws.send(JSON.stringify({ action: 'identified', parameters: { 'predefs': refList } }));

        break;
      case "useRef":
        const refName = parameters.refName;
        ws.reffa = `/opt/shiba2motif/predefinedref/${refName}`;
        console.log(ws.reffa);
        ws.send(JSON.stringify({ action: 'predefinedRef', parameters: { 'info': refName } }));
        break;


    }

  });

});


http.createServer(function (req, res) {


  const ip = req.socket.remoteAddress;
  let ws = null;
  const wssC = Array.from(wss.clients);
  // for the request, find the ws that has the same ip as the request
  for (const wsClient in wssC) {
    console.log(ip);
    console.log(wssC[wsClient]._socket.remoteAddress);
    if (wssC[wsClient]._socket.remoteAddress == ip) {
      ws = wssC[wsClient];
    }
  }

  console.log('received file from ' + ip);
  if (req.url == '/five') {
    var form = new formidable.IncomingForm();
    form.parse(req, function (err, fields, files) {
      var oldpath = files.filename.filepath;
      var newpath = `/tmp/input${ip}/five.txt`;
      fs.rename(oldpath, newpath, function (err) {
        if (err) throw err;
        res.write('File uploaded and moved!');
        res.end();
        ws.send(JSON.stringify({ action: 'filesystemChange', parameters: { 'info': newpath } }));
      });
    });
  }

  if (req.url == '/mxe') {
    var form = new formidable.IncomingForm();
    form.parse(req, function (err, fields, files) {
      var oldpath = files.filename.filepath;
      var newpath = `/tmp/input${ip}/mxe.txt`;
      fs.rename(oldpath, newpath, function (err) {
        if (err) throw err;
        res.write('File uploaded and moved!');
        res.end();
        ws.send(JSON.stringify({ action: 'filesystemChange', parameters: { 'info': newpath } }));
      });
    });
  }

  if (req.url == '/ri') {
    var form = new formidable.IncomingForm();
    form.parse(req, function (err, fields, files) {
      var oldpath = files.filename.filepath;
      var newpath = `/tmp/input${ip}/ri.txt`;
      fs.rename(oldpath, newpath, function (err) {
        if (err) throw err;
        res.write('File uploaded and moved!');
        res.end();
        ws.send(JSON.stringify({ action: 'filesystemChange', parameters: { 'info': newpath } }));
      });
    });
  }

  if (req.url == '/se') {
    var form = new formidable.IncomingForm();
    form.parse(req, function (err, fields, files) {
      var oldpath = files.filename.filepath;
      var newpath = `/tmp/input${ip}/se.txt`;
      fs.rename(oldpath, newpath, function (err) {
        if (err) throw err;
        res.write('File uploaded and moved!');
        res.end();
        ws.send(JSON.stringify({ action: 'filesystemChange', parameters: { 'info': newpath } }));
      });
    });
  }

  if (req.url == '/three') {
    var form = new formidable.IncomingForm();
    form.parse(req, function (err, fields, files) {
      var oldpath = files.filename.filepath;
      var newpath = `/tmp/input${ip}/three.txt`;
      fs.rename(oldpath, newpath, function (err) {
        if (err) throw err;
        res.write('File uploaded and moved!');
        res.end();
        ws.send(JSON.stringify({ action: 'filesystemChange', parameters: { 'info': newpath } }));
      });
    });
  }

  if (req.url == '/ref') {
    var form = new formidable.IncomingForm();
    form.parse(req, function (err, fields, files) {
      console.log(files)
      var oldpath = files.filename.filepath;
      var newpath = `/tmp/ref${ip}/ref.fa`;
      fs.rename(oldpath, newpath, function (err) {
        if (err) throw err;
        res.write('File uploaded and moved!');
        res.end();
        ws.send(JSON.stringify({ action: 'filesystemChange', parameters: { 'info': newpath } }));
      });
    });
  }
}).listen(8080);
