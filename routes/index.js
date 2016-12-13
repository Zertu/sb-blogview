const express = require('express'),
  router = express.Router(),
  user = require('./user'),
  session = require('express-session')
  /* GET home page. */
let app= express()
app.use(session({
    secret: 'mysbblog', // 建议使用 128 个字符的随机字符串
    cookie: { maxAge: 1000 * 60 * 60 * 24 }
}))
router.get('/', function (req, res, next) {
  res.render('index', {
    title: 'sb-blog'
  });
});
router.get('/article',(req, res) => {
  res.render('article', {
    title: 'sb-blog',session:req.session
  })
})


router.post('/login', (req, res) => {
  user.checklogin(req, res)
})

module.exports = router;