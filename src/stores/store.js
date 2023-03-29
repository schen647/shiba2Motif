import { ref, computed } from 'vue'
import { defineStore } from 'pinia'

const log1 = ref('')
const availableRef =ref([])

// establish a ws to 192.168.1.114:8080
const socket = new WebSocket("ws://10.10.10.167:8089");
socket.onopen = function (e) {
  sendWS({
    action: 'ident',
    parameters: { role: 'client' }
  })
};

socket.onmessage = function (event) {
  //netResp.value = event.data.replace(/\n/g, '</br>');
  // netResp.value = event.data
  
  const data = JSON.parse(event.data)
  const action = data.action
  const parameters = data.parameters
  console.log(data)
  switch (action) {
    case 'identified':
      for (const ref1 in parameters.predefs) {
        availableRef.value.push(parameters.predefs[ref1])
      }
      break;
    
    case 'filesystemChange':
      log1.value += JSON.stringify(data)+'\n'
      break;

    case 'predefinedRef':
      log1.value += JSON.stringify(data)+'\n'
      break;


  }
};

function runInfer() {
  sendWS({
    action: 'run',
    parameters: {}
  })
}

function useRef(refName) {
  sendWS({
    action: 'useRef',
    parameters: {refName}
  })
}

socket.onclose = function (event) {
  if (event.wasClean) {
    alert(`[close] Connection closed cleanly, code=${event.code} reason=${event.reason}`);
  } else {
    // e.g. server process killed or network down
    // event.code is usually 1006 in this case
    alert('[close] Connection died');
  }
};

socket.onerror = function (error) {
  alert(`[error]`);
};

function sendWS(msg) {
  console.log('sending', JSON.stringify(msg))
  socket.send(JSON.stringify(msg));
}


export const useUserStore = defineStore('user', () => {
  return { runInfer, log1, availableRef, useRef }
})
