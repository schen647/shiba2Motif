<script>
import { RouterLink, RouterView } from 'vue-router'
import { mapActions, mapState } from 'pinia'
import { useUserStore } from './stores/store'

export default {
  // feed those
  // props: ['channels', 'chatLog'],
  data() {
    return {
      username: '',
      userName2beSubmitted: '',
      password: '',
      activeSelector: 'ass',
    }
  },

  computed: {
    //...mapState(useUserStore, ['netRespAssist','loggedIn']),

  },
  mounted() {

  },
  updated() {
  },
  methods: {
    //...mapActions(useUserStore, [ 'login']),
    setActiveSelector(se) {
      this.activeSelector = se;
    },
    submitUsrName() {
      this.userName2beSubmitted = this.username

      console.log('received username')
    },
    clearunsuccessfulLogin() {
      if (this.loggedIn)
        return
      this.userName2beSubmitted = ''
    },
    logMeIn() {
      this.login(this.userName2beSubmitted, this.password)
    }

  },
}
</script>

<template>
  <header>

  </header>



  <RouterLink to="/findMotif">
    findMotif
  </RouterLink>

  <RouterView />
</template>

<style scoped></style>
