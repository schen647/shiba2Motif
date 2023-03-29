import { createRouter, createWebHistory } from 'vue-router'
import findMotif from '../views/findMotif.vue'

const router = createRouter({
  history: createWebHistory(import.meta.env.BASE_URL),
  routes: [

    {
      path: '/findMotif',
      name: 'findMotif',
      component: findMotif
    },
   
  ]
})

export default router
