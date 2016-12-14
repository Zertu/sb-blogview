const express = require('express'),
  router = express.Router(),
  user = require('./user')
  /* GET home page. */
router.get('/', function (req, res, next) {
  res.render('index', {
    title: 'sb-blog'
  });
})

router.get('/article', (req, res) => {
  res.render('article', {
    title: 'sb-blog',
    session: req.session.user
  })
})

router.get('/logout',(req,res)=>{
  req.session.user=null
  res.render('article',{
    title:'sb-blog',
    session:req.session.user
  })
})

// post 请求处理
router.post('/login', (req, res) => {
  if(req.session.user){
    res.render('article', {
      title: 'sb-blog',
      session: req.session.user
    })
  }
  let row = user.checklogin(req, res, function (rows) {
    res.render('article', {
      title: 'sb-blog',
      session: req.session.user
    })
  })
})

module.exports = router;