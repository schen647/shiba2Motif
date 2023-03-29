<script>
import { mapActions, mapState } from 'pinia'
import { useUserStore } from '../stores/store'
export default {
  // feed those
  // props: ['channels', 'chatLog'],
  data() {
    return {
      useCustomRef: false,

    }
  },

  computed: {
    ...mapState(useUserStore, ['log1', 'availableRef']),
    netProcess() {
      function removeLastn(input) {
        var str = input
        var word = '\n';
        var newWord = '';

        // find the index of last time word was used
        // please note lastIndexOf() is case sensitive
        var n = str.lastIndexOf(word);

        // slice the string in 2, one from the start to the lastIndexOf
        // and then replace the word in the rest
        str = str.slice(0, n) + str.slice(n).replace(word, newWord);
        return str
      }
      // split this.userinput into 2 parts: the last 1024 letters and all the previous letters

      return this.userInput + this.netRespAssist
    }
  },
  mounted() {
    //let currentInput = ''


  },
  unmounted() {
    //clearInterval(this.querry)
  },
  updated() {
  },
  methods: {
    ...mapActions(useUserStore, ['runInfer', 'useRef']),
    run() {
      this.runInfer()
    },
    switchRef() {
      this.useCustomRef = !this.useCustomRef
    },
    useThisRef(nam){
      this.useRef(nam)
    }


  },
}
</script>
<template>
  find motif
  <pre style="color:red;right:0;position:absolute;">{{ String(log1) }}</pre>
  <iframe name="dummyframe" id="dummyframe" style="display: none;"></iframe>
  <div @click="switchRef">switch ref</div>
  <form v-if="useCustomRef" action="http://10.10.10.167:8080/ref" target="dummyframe" enctype="multipart/form-data" method="post">
    <p>ref</p>
    <input type="file" id="myFile" name="filename">
    <input type="submit">

  </form>
  <div v-else>
    <div v-for="ref1 in availableRef">
    <div @click="useThisRef(ref1)">
      {{ ref1 }}
    </div>
    
    </div>

  </div>

  <form action="http://10.10.10.167:8080/three" target="dummyframe" enctype="multipart/form-data" method="post">
    <p>three</p>
    <input type="file" id="myFile" name="filename">
    <input type="submit">

  </form>

  <form action="http://10.10.10.167:8080/se" target="dummyframe" enctype="multipart/form-data" method="post">
    <p>se</p>
    <input type="file" id="myFile" name="filename">
    <input type="submit">

  </form>

  <form action="http://10.10.10.167:8080/ri" target="dummyframe" enctype="multipart/form-data" method="post">
    <p>ri</p>
    <input type="file" id="myFile" name="filename">
    <input type="submit">

  </form>

  <form action="http://10.10.10.167:8080/mxe" target="dummyframe" enctype="multipart/form-data" method="post">
    <p>mxe</p>
    <input type="file" id="myFile" name="filename">
    <input type="submit">

  </form>

  <form action="http://10.10.10.167:8080/five" target="dummyframe" enctype="multipart/form-data" method="post">
    <p>five</p>
    <input type="file" id="myFile" name="filename">
    <input type="submit">

  </form>
  <div @click="run()">RUN</div>
</template>

<style></style>
