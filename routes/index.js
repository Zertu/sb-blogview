const express = require('express')
const router = express.Router()
const user = require('./user')
/* GET home page. */
router.get('/', function(req, res, next) {
  res.render('index', { title: 'sb-blog' });
});
router.get('/article',(req,res)=>{
  res.render('article',{ title: 'sb-blog' })
})


router.post('/login',(req,res)=>{
  user.checklogin(req,res)
})

module.exports = router;
